# =============================================================================
# binary_response_proportions.R
# Dataset: b11  |  Colonne: c2.2. … c2.27.
# Categorie risposta: 1, 0, NC, NS
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

# ── PARAMETRI ─────────────────────────────────────────────────────────────────
p0 <- 0.50   # ipotesi nulla per i test (modificare se necessario)

# ── IDENTIFICAZIONE COLONNE c2.X. ─────────────────────────────────────────────
colonne_c2 <- grep("^c2\\.", names(b11), value = TRUE)
cat(sprintf("Colonne trovate (%d): %s\n\n",
            length(colonne_c2), paste(colonne_c2, collapse = ", ")))

# ── LOOP SU OGNI COLONNA ──────────────────────────────────────────────────────
risultati_proporzioni <- list()
risultati_test        <- list()

for (col in colonne_c2) {

  risposte <- as.character(b11[[col]])
  N_tot    <- length(risposte)

  # Proporzioni 4 categorie
  freq <- table(factor(risposte, levels = c("1", "0", "NC", "NS")))
  prop <- prop.table(freq)

  risultati_proporzioni[[col]] <- data.frame(
    colonna  = col,
    N_totale = N_tot,
    n_1      = as.integer(freq["1"]),
    n_0      = as.integer(freq["0"]),
    n_NC     = as.integer(freq["NC"]),
    n_NS     = as.integer(freq["NS"]),
    pct_1    = round(prop["1"]  * 100, 1),
    pct_0    = round(prop["0"]  * 100, 1),
    pct_NC   = round(prop["NC"] * 100, 1),
    pct_NS   = round(prop["NS"] * 100, 1),
    row.names = NULL
  )

  # Sottoinsieme valido (escludi NC e NS)
  validi    <- risposte[risposte %in% c("1", "0")]
  n_val     <- length(validi)
  n_pos     <- sum(validi == "1")
  n_esclusi <- N_tot - n_val

  if (n_val < 5) {
    risultati_test[[col]] <- data.frame(
      colonna       = col, N_validi = n_val, n_esclusi = n_esclusi,
      n_1 = n_pos,  n_0 = n_val - n_pos,
      prop_1 = NA, ic95_lower = NA, ic95_upper = NA,
      p_binom = NA, p_prop = NA, significativo = NA,
      row.names = NULL
    )
    next
  }

  bt <- binom.test(n_pos, n_val, p = p0, alternative = "two.sided", conf.level = 0.95)
  pt <- prop.test(n_pos, n_val,  p = p0, alternative = "two.sided",
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
    p_binom       = round(bt$p.value, 4),
    p_prop        = round(pt$p.value, 4),
    significativo = ifelse(bt$p.value < 0.05, "SI *", "no"),
    row.names = NULL
  )
}

# ── TABELLE RIASSUNTIVE ───────────────────────────────────────────────────────
tab_prop <- bind_rows(risultati_proporzioni)
tab_test <- bind_rows(risultati_test)

cat("\n============================================================\n")
cat("  PROPORZIONI – TUTTE E 4 LE CATEGORIE\n")
cat("============================================================\n\n")
print(tab_prop, row.names = FALSE)

cat("\n============================================================\n")
cat(sprintf("  TEST (H0: p = %.2f) – solo rispondenti validi 0/1\n", p0))
cat("============================================================\n\n")
print(tab_test, row.names = FALSE)

sig <- tab_test$colonna[!is.na(tab_test$significativo) & tab_test$significativo == "SI *"]
cat(sprintf("\nColonne significative (p < 0.05): %s\n",
            if (length(sig) == 0) "nessuna" else paste(sig, collapse = ", ")))

# ── EXPORT (decommentare se necessario) ───────────────────────────────────────
# write.csv(tab_prop, "proporzioni_c2.csv", row.names = FALSE)
# write.csv(tab_test, "test_binari_c2.csv", row.names = FALSE)

# ── GRAFICI ───────────────────────────────────────────────────────────────────

# 1. Barre impilate – distribuzione 4 categorie per colonna
tab_prop %>%
  select(colonna, pct_1, pct_0, pct_NC, pct_NS) %>%
  pivot_longer(-colonna, names_to = "cat", values_to = "pct") %>%
  mutate(
    cat     = factor(cat,
                     levels = c("pct_1", "pct_0", "pct_NC", "pct_NS"),
                     labels = c("1 (Sì)", "0 (No)", "NC", "NS")),
    colonna = factor(colonna, levels = colonne_c2)
  ) %>%
  ggplot(aes(x = colonna, y = pct, fill = cat)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c("1 (Sì)" = "#2166ac", "0 (No)" = "#d73027",
                                "NC"     = "#a6d96a", "NS"     = "#fdae61")) +
  scale_y_continuous(labels = label_percent(scale = 1)) +
  labs(title = "Distribuzione risposte per colonna",
       x = NULL, y = "%", fill = NULL) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")

# 2. Proporzione "1" con IC 95% (solo validi)
tab_test %>%
  filter(!is.na(prop_1)) %>%
  mutate(colonna = factor(colonna, levels = colonne_c2)) %>%
  ggplot(aes(x = colonna, y = prop_1, colour = significativo)) +
  geom_hline(yintercept = p0, linetype = "dashed", colour = "grey50") +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = ic95_lower, ymax = ic95_upper), width = 0.25) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_colour_manual(values = c("SI *" = "#d73027", "no" = "#4393c3"),
                      name = "p < 0.05") +
  labs(title = "Proporzione risposte '1' per colonna (IC 95%)",
       subtitle = sprintf("Linea tratteggiata = H0: p = %.2f  |  Esclusi NC e NS", p0),
       x = NULL, y = "Proporzione '1'") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
