# =============================================================================
# binary_response_proportions.R
# -----------------------------------------------------------------------------
# Analisi di risposte con 4 categorie: 1, 0, NC, NS
# Colonne: c2.2. c2.3. c2.5. c2.7. c2.9. c2.11. c2.13. c2.15.
#          c2.17. c2.19. c2.21. c2.23. c2.25. c2.27.
#
# Flusso per ogni colonna:
#   1. Proporzioni su tutte e 4 le categorie
#   2. Esclusione NC e NS
#   3. binom.test + prop.test (H0: p = 0.50) sul sottoinsieme 0/1
#   4. Tabella riassuntiva unica + grafici
# =============================================================================


# ── 0. PACCHETTI ──────────────────────────────────────────────────────────────

# install.packages(c("dplyr", "tidyr", "ggplot2", "scales"))
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)


# ── 1. CARICAMENTO DATI ───────────────────────────────────────────────────────
# Decommentare la riga adatta al proprio formato:

# dati <- read.csv("file.csv", stringsAsFactors = FALSE)
# dati <- readxl::read_excel("file.xlsx")

# ── DATI DI ESEMPIO (rimuovere quando si caricano i dati reali) ───────────────
set.seed(42)
n <- 200
cols_c2 <- c("c2.2.", "c2.3.", "c2.5.", "c2.7.", "c2.9.", "c2.11.",
             "c2.13.", "c2.15.", "c2.17.", "c2.19.", "c2.21.", "c2.23.",
             "c2.25.", "c2.27.")

dati <- as.data.frame(
  lapply(cols_c2, function(x)
    sample(c("1", "0", "NC", "NS"), n, replace = TRUE,
           prob = c(0.45, 0.30, 0.15, 0.10))
  )
)
names(dati) <- cols_c2
# ─────────────────────────────────────────────────────────────────────────────


# ── 2. IDENTIFICAZIONE AUTOMATICA DELLE COLONNE c2.X. ────────────────────────

colonne_c2 <- grep("^c2\\.", names(dati), value = TRUE)
cat(sprintf("Colonne analizzate (%d): %s\n\n",
            length(colonne_c2), paste(colonne_c2, collapse = ", ")))

# Ipotesi nulla per il test (modificare se necessario)
p0 <- 0.50


# ── 3. LOOP SU OGNI COLONNA ───────────────────────────────────────────────────

# Contenitori per i risultati
risultati_proporzioni <- list()   # proporzioni 4 categorie
risultati_test        <- list()   # risultati test binari

for (col in colonne_c2) {

  risposte <- dati[[col]]
  N_tot    <- length(risposte)

  # ── 3a. Proporzioni 4 categorie ──
  freq <- table(factor(risposte, levels = c("1", "0", "NC", "NS")))
  prop <- prop.table(freq)

  risultati_proporzioni[[col]] <- data.frame(
    colonna     = col,
    N_totale    = N_tot,
    n_1         = as.integer(freq["1"]),
    n_0         = as.integer(freq["0"]),
    n_NC        = as.integer(freq["NC"]),
    n_NS        = as.integer(freq["NS"]),
    pct_1       = round(prop["1"]  * 100, 1),
    pct_0       = round(prop["0"]  * 100, 1),
    pct_NC      = round(prop["NC"] * 100, 1),
    pct_NS      = round(prop["NS"] * 100, 1)
  )

  # ── 3b. Sottoinsieme valido (0/1) ──
  validi  <- risposte[risposte %in% c("1", "0")]
  n_val   <- length(validi)
  n_pos   <- sum(validi == "1")
  n_esclusi <- N_tot - n_val

  if (n_val < 5) {
    # Troppo poche osservazioni: segnala e salta il test
    risultati_test[[col]] <- data.frame(
      colonna     = col,
      N_validi    = n_val,
      n_esclusi   = n_esclusi,
      n_1         = n_pos,
      n_0         = n_val - n_pos,
      prop_1      = NA_real_,
      ic95_lower  = NA_real_,
      ic95_upper  = NA_real_,
      p_value_binom = NA_real_,
      p_value_prop  = NA_real_,
      significativo = NA
    )
    next
  }

  # binom.test
  bt <- binom.test(n_pos, n_val, p = p0, alternative = "two.sided", conf.level = 0.95)
  # prop.test
  pt <- prop.test(n_pos, n_val, p = p0, alternative = "two.sided",
                  conf.level = 0.95, correct = TRUE)

  risultati_test[[col]] <- data.frame(
    colonna       = col,
    N_validi      = n_val,
    n_esclusi     = n_esclusi,
    n_1           = n_pos,
    n_0           = n_val - n_pos,
    prop_1        = round(n_pos / n_val, 4),
    ic95_lower    = round(bt$conf.int[1], 4),
    ic95_upper    = round(bt$conf.int[2], 4),
    p_value_binom = round(bt$p.value, 4),
    p_value_prop  = round(pt$p.value, 4),
    significativo = ifelse(bt$p.value < 0.05, "SI *", "no")
  )
}


# ── 4. TABELLE RIASSUNTIVE ────────────────────────────────────────────────────

tab_prop <- bind_rows(risultati_proporzioni)
tab_test <- bind_rows(risultati_test)

cat("\n============================================================\n")
cat("  PROPORZIONI – TUTTE E 4 LE CATEGORIE\n")
cat("============================================================\n\n")
print(tab_prop, row.names = FALSE, digits = 3)

cat("\n============================================================\n")
cat(sprintf("  TEST BINOMIALE / PROP.TEST  (H0: p = %.2f)\n", p0))
cat("  Solo rispondenti validi (esclusi NC e NS)\n")
cat("============================================================\n\n")
print(tab_test, row.names = FALSE, digits = 4)

# Riepilogo colonne significative
sig <- tab_test[!is.na(tab_test$significativo) & tab_test$significativo == "SI *", "colonna"]
cat(sprintf("\nColonne con p < 0.05 (%d): %s\n",
            length(sig),
            if (length(sig) == 0) "nessuna" else paste(sig, collapse = ", ")))


# ── 5. EXPORT RISULTATI ───────────────────────────────────────────────────────
# write.csv(tab_prop, "proporzioni_c2.csv",  row.names = FALSE)
# write.csv(tab_test, "test_binari_c2.csv",  row.names = FALSE)


# ── 6. GRAFICI ────────────────────────────────────────────────────────────────

# ── 6a. Heatmap proporzioni 4 categorie ──────────────────────────────────────
plot_prop <- tab_prop %>%
  select(colonna, pct_1, pct_0, pct_NC, pct_NS) %>%
  pivot_longer(-colonna, names_to = "categoria", values_to = "pct") %>%
  mutate(
    categoria = factor(categoria,
                       levels = c("pct_1", "pct_0", "pct_NC", "pct_NS"),
                       labels = c("1 (Sì)", "0 (No)", "NC", "NS")),
    colonna   = factor(colonna, levels = rev(colonne_c2))
  )

ggplot(plot_prop, aes(x = categoria, y = colonna, fill = pct)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.1f%%", pct)), size = 3) +
  scale_fill_gradient(low = "#f7fbff", high = "#2171b5",
                      name = "%", limits = c(0, 100)) +
  labs(title = "Distribuzione risposte per colonna (%)",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(axis.text.y = element_text(size = 10),
        panel.grid  = element_blank())


# ── 6b. Proporzione risposte '1' (validi) con IC 95% ─────────────────────────
plot_test <- tab_test %>%
  filter(!is.na(prop_1)) %>%
  mutate(colonna = factor(colonna, levels = colonne_c2))

ggplot(plot_test, aes(x = colonna, y = prop_1, colour = significativo)) +
  geom_hline(yintercept = p0, linetype = "dashed", colour = "grey50") +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = ic95_lower, ymax = ic95_upper), width = 0.25) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_colour_manual(values = c("SI *" = "#d73027", "no" = "#4393c3"),
                      name = sprintf("p < 0.05\n(H0: p=%.2f)", p0)) +
  labs(title = "Proporzione risposte '1' per colonna (IC 95%)",
       subtitle = "Punti rossi = significativo; linea tratteggiata = ipotesi nulla",
       x = "Colonna", y = "Proporzione '1'") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# ── 6c. Barre impilate 4 categorie ───────────────────────────────────────────
plot_prop_stacked <- plot_prop %>%
  mutate(colonna = factor(colonna, levels = colonne_c2))

ggplot(plot_prop_stacked, aes(x = colonna, y = pct, fill = categoria)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c("1 (Sì)" = "#2166ac", "0 (No)" = "#d73027",
                                "NC"     = "#a6d96a", "NS"     = "#fdae61")) +
  scale_y_continuous(labels = label_percent(scale = 1)) +
  labs(title = "Distribuzione risposte per colonna (barre impilate)",
       x = "Colonna", y = "%", fill = NULL) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")

cat("\nScript completato.\n")
