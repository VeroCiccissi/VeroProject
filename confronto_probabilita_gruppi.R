# =============================================================================
# confronto_probabilita_gruppi.R
# Dataset vaccinate: b11   |  Dataset non vaccinate: SubNV
# Variabili: g1, g2, g3, g4  (scala 0-100, NC = non risposta)
#
# Flusso:
#   1. Rimozione NC e conversione numerica
#   2. Statistiche descrittive per gruppo
#   3. Test di Wilcoxon/Mann-Whitney (confronto distribuzioni)
#   4. Grafici boxplot + densità
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)

variabili <- c("g1", "g2", "g3", "g4")

# ── 1. PULIZIA: rimuovi NC e converti a numerico ──────────────────────────────

pulisci <- function(df, vars, gruppo) {
  df %>%
    select(all_of(vars)) %>%
    mutate(across(everything(), ~ as.numeric(ifelse(. == "NC", NA, as.character(.))))) %>%
    mutate(gruppo = gruppo)
}

vacc   <- pulisci(b11,   variabili, "Vaccinate")
nonvac <- pulisci(SubNV, variabili, "Non vaccinate")

combined <- bind_rows(vacc, nonvac) %>%
  mutate(gruppo = factor(gruppo, levels = c("Vaccinate", "Non vaccinate")))

cat(sprintf("Vaccinate:     %d osservazioni\n", nrow(vacc)))
cat(sprintf("Non vaccinate: %d osservazioni\n", nrow(nonvac)))

# NC esclusi per variabile
cat("\nNC esclusi (non entrano nell'analisi):\n")
for (v in variabili) {
  nc_v  <- sum(b11[[v]]   == "NC", na.rm = TRUE)
  nc_nv <- sum(SubNV[[v]] == "NC", na.rm = TRUE)
  cat(sprintf("  %s → vaccinate: %d NC  |  non vaccinate: %d NC\n", v, nc_v, nc_nv))
}


# ── 2. STATISTICHE DESCRITTIVE ────────────────────────────────────────────────

cat("\n============================================================\n")
cat("  STATISTICHE DESCRITTIVE (esclusi NC)\n")
cat("============================================================\n\n")

desc <- combined %>%
  pivot_longer(all_of(variabili), names_to = "variabile", values_to = "valore") %>%
  filter(!is.na(valore)) %>%
  group_by(variabile, gruppo) %>%
  summarise(
    n      = n(),
    media  = round(mean(valore), 1),
    mediana= round(median(valore), 1),
    IQR    = round(IQR(valore), 1),
    sd     = round(sd(valore), 1),
    min    = min(valore),
    max    = max(valore),
    .groups = "drop"
  ) %>%
  arrange(variabile, gruppo)

print(desc, row.names = FALSE)


# ── 3. TEST DI WILCOXON / MANN-WHITNEY ───────────────────────────────────────

cat("\n============================================================\n")
cat("  WILCOXON / MANN-WHITNEY  (H0: distribuzioni uguali)\n")
cat("============================================================\n\n")

risultati_wilcox <- list()

for (v in variabili) {

  x_v  <- vacc[[v]][!is.na(vacc[[v]])]
  x_nv <- nonvac[[v]][!is.na(nonvac[[v]])]

  wt <- wilcox.test(x_v, x_nv,
                    alternative = "two.sided",
                    conf.int    = TRUE,
                    conf.level  = 0.95)

  risultati_wilcox[[v]] <- data.frame(
    variabile        = v,
    n_vacc           = length(x_v),
    n_nonvacc        = length(x_nv),
    mediana_vacc     = round(median(x_v), 1),
    mediana_nonvacc  = round(median(x_nv), 1),
    diff_mediane     = round(median(x_v) - median(x_nv), 1),
    stima_HL         = round(wt$estimate, 2),   # Hodges-Lehmann
    ic95_lower       = round(wt$conf.int[1], 2),
    ic95_upper       = round(wt$conf.int[2], 2),
    p_value          = round(wt$p.value, 4),
    significativo    = ifelse(wt$p.value < 0.05, "SI *", "no"),
    row.names = NULL
  )
}

tab_wilcox <- bind_rows(risultati_wilcox)
print(tab_wilcox, row.names = FALSE)

sig_w <- tab_wilcox$variabile[tab_wilcox$significativo == "SI *"]
cat(sprintf("\nVariabili significative (p < 0.05): %s\n",
            if (length(sig_w) == 0) "nessuna" else paste(sig_w, collapse = ", ")))

# ── EXPORT (decommentare se necessario) ───────────────────────────────────────
# write.csv(desc,       "descrittive_g1g4.csv",  row.names = FALSE)
# write.csv(tab_wilcox, "wilcoxon_g1g4.csv",     row.names = FALSE)


# ── 4. GRAFICI ────────────────────────────────────────────────────────────────

dati_long <- combined %>%
  pivot_longer(all_of(variabili), names_to = "variabile", values_to = "valore") %>%
  filter(!is.na(valore))

# Etichette n per ogni pannello
etichette_n <- dati_long %>%
  group_by(variabile, gruppo) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(label = paste0(gruppo, "\n(n=", n, ")"))

# ── Grafico 1: Boxplot affiancati per variabile ───────────────────────────────
ggplot(dati_long, aes(x = gruppo, y = valore, fill = gruppo)) +
  geom_boxplot(outlier.size = 1, alpha = 0.7, width = 0.5) +
  geom_jitter(width = 0.15, alpha = 0.15, size = 0.8) +
  facet_wrap(~ variabile, nrow = 1) +
  scale_fill_manual(values = c("Vaccinate" = "#2166ac", "Non vaccinate" = "#d73027")) +
  geom_text(data = tab_wilcox,
            aes(x = 1.5, y = Inf,
                label = sprintf("p = %.4f%s", p_value,
                                ifelse(significativo == "SI *", " *", ""))),
            inherit.aes = FALSE, vjust = 1.5, size = 3.5) +
  labs(title = "Distribuzione valori per variabile e gruppo",
       subtitle = "Esclusi NC  |  * p < 0.05 (Wilcoxon)",
       x = NULL, y = "Valore (0–100)", fill = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 20, hjust = 1))

# ── Grafico 2: Densità sovrapposte ────────────────────────────────────────────
ggplot(dati_long, aes(x = valore, fill = gruppo, colour = gruppo)) +
  geom_density(alpha = 0.35, linewidth = 0.8) +
  facet_wrap(~ variabile, nrow = 2, scales = "free_y") +
  scale_fill_manual(values   = c("Vaccinate" = "#2166ac", "Non vaccinate" = "#d73027")) +
  scale_colour_manual(values = c("Vaccinate" = "#2166ac", "Non vaccinate" = "#d73027")) +
  labs(title = "Distribuzione densità per variabile e gruppo",
       subtitle = "Esclusi NC",
       x = "Valore (0–100)", y = "Densità", fill = NULL, colour = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")

cat("\nScript completato.\n")
