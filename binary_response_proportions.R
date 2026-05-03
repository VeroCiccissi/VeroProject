# =============================================================================
# binary_response_proportions.R
# -----------------------------------------------------------------------------
# Analisi di risposte con 4 categorie: 1, 0, NC (Non Classificabile), NS (Non Sa)
#
# Flusso:
#   1. Dati di esempio con 4 categorie di risposta
#   2. Proporzioni su tutte e 4 le categorie (incluse NC/NS)
#   3. Esclusione di NC e NS  →  analisi solo sui rispondenti validi (0/1)
#   4. Test statistici su risposte binarie: prop.test, binom.test, chi-quadro
# =============================================================================


# ── 0. PACCHETTI ──────────────────────────────────────────────────────────────

# install.packages(c("dplyr", "ggplot2", "scales"))   # eseguire solo se mancanti
library(dplyr)
library(ggplot2)
library(scales)


# ── 1. DATI DI ESEMPIO ────────────────────────────────────────────────────────
# Sostituire 'risposte' con il proprio vettore / colonna del data frame

set.seed(42)
n <- 200

dati <- data.frame(
  id       = 1:n,
  gruppo   = sample(c("A", "B"), n, replace = TRUE),
  risposta = sample(
    c("1", "0", "NC", "NS"),
    n,
    replace = TRUE,
    prob    = c(0.45, 0.30, 0.15, 0.10)   # probabilità attese
  )
)

# Per caricare dati reali usare, per esempio:
# dati <- read.csv("file.csv")
# oppure: dati <- readxl::read_excel("file.xlsx")


# ── 2. PROPORZIONI – TUTTE E 4 LE CATEGORIE ──────────────────────────────────

cat("\n====================================================\n")
cat("  STEP 1: Distribuzione completa (4 categorie)\n")
cat("====================================================\n\n")

tab_completa <- dati %>%
  count(risposta, name = "n") %>%
  mutate(
    proporzione = n / sum(n),
    percentuale = proporzione * 100,
    categoria   = case_when(
      risposta == "1"  ~ "Sì / Presente",
      risposta == "0"  ~ "No / Assente",
      risposta == "NC" ~ "Non Classificabile",
      risposta == "NS" ~ "Non Sa / Non Risponde",
      TRUE             ~ risposta
    )
  ) %>%
  arrange(match(risposta, c("1", "0", "NC", "NS")))

print(tab_completa, digits = 3)

cat(sprintf(
  "\nTotale osservazioni: %d\n  di cui NC:  %d (%.1f%%)\n  di cui NS:  %d (%.1f%%)\n",
  nrow(dati),
  sum(dati$risposta == "NC"), mean(dati$risposta == "NC") * 100,
  sum(dati$risposta == "NS"), mean(dati$risposta == "NS") * 100
))


# ── 3. ESCLUSIONE NC / NS ─────────────────────────────────────────────────────

cat("\n====================================================\n")
cat("  STEP 2: Esclusione NC e NS → analisi su 0 e 1\n")
cat("====================================================\n\n")

dati_binari <- dati %>%
  filter(risposta %in% c("1", "0")) %>%
  mutate(risposta_num = as.integer(risposta))   # 1 → 1, 0 → 0

n_esclusi <- nrow(dati) - nrow(dati_binari)
cat(sprintf("Osservazioni escluse (NC + NS): %d\n", n_esclusi))
cat(sprintf("Osservazioni valide per l'analisi: %d\n\n", nrow(dati_binari)))

# Conteggi base
n_val  <- nrow(dati_binari)
n_pos  <- sum(dati_binari$risposta_num == 1)
n_neg  <- sum(dati_binari$risposta_num == 0)
p_obs  <- n_pos / n_val

cat(sprintf("  Risposte '1':  %d  (%.1f%%)\n", n_pos, p_obs * 100))
cat(sprintf("  Risposte '0':  %d  (%.1f%%)\n", n_neg, (1 - p_obs) * 100))


# ── 4. TEST STATISTICI ────────────────────────────────────────────────────────

cat("\n====================================================\n")
cat("  STEP 3: Test statistici sulle risposte binarie\n")
cat("====================================================\n\n")


# 4a. binom.test  – proporzione osservata vs ipotesi nulla p0 = 0.50
# (modificare p0 secondo l'ipotesi nulla di interesse)
p0 <- 0.50

cat("--- Test binomiale esatto  (H0: p =", p0, ") ---\n")
bt <- binom.test(n_pos, n_val, p = p0, alternative = "two.sided", conf.level = 0.95)
print(bt)


# 4b. prop.test  – alternativa normale (utile con campioni grandi)
cat("\n--- prop.test  (H0: p =", p0, ") ---\n")
pt <- prop.test(n_pos, n_val, p = p0, alternative = "two.sided",
                conf.level = 0.95, correct = TRUE)
print(pt)

cat(sprintf(
  "\nProporzione stimata:  %.4f  (IC 95%%: %.4f – %.4f)\n",
  p_obs, pt$conf.int[1], pt$conf.int[2]
))


# 4c. Confronto proporzionale tra due gruppi (A vs B)
cat("\n--- Confronto tra gruppi (A vs B)  – prop.test  ---\n")

tab_gruppi <- dati_binari %>%
  group_by(gruppo) %>%
  summarise(
    n_tot  = n(),
    n_pos  = sum(risposta_num == 1),
    prop   = n_pos / n_tot,
    .groups = "drop"
  )
print(tab_gruppi)

pt2 <- prop.test(
  x       = tab_gruppi$n_pos,
  n       = tab_gruppi$n_tot,
  alternative = "two.sided",
  conf.level  = 0.95,
  correct     = TRUE
)
print(pt2)

cat(sprintf(
  "\nDifferenza di proporzione (A − B): %.4f  (IC 95%%: %.4f – %.4f)\n",
  diff(rev(pt2$estimate)),    # A - B
  pt2$conf.int[1],
  pt2$conf.int[2]
))


# 4d. Test chi-quadro su tabella di contingenza (include verifica assunzioni)
cat("\n--- Chi-quadro su tabella 2×2 (risposta × gruppo) ---\n")

tab_2x2 <- table(dati_binari$risposta, dati_binari$gruppo)
print(tab_2x2)

# Verifica celle attese >= 5
expected <- chisq.test(tab_2x2)$expected
if (any(expected < 5)) {
  cat("ATTENZIONE: celle attese < 5, usare test di Fisher.\n")
  fisher <- fisher.test(tab_2x2, conf.level = 0.95)
  print(fisher)
} else {
  chi2 <- chisq.test(tab_2x2, correct = FALSE)
  print(chi2)
  cat(sprintf("OR (Fisher): %.3f\n", fisher.test(tab_2x2)$estimate))
}


# ── 5. GRAFICO ────────────────────────────────────────────────────────────────

cat("\n====================================================\n")
cat("  STEP 4: Grafici\n")
cat("====================================================\n\n")

# 5a. Grafico a barre – distribuzione completa con 4 categorie
tab_plot <- tab_completa %>%
  mutate(risposta = factor(risposta, levels = c("1", "0", "NC", "NS")))

p1 <- ggplot(tab_plot, aes(x = risposta, y = percentuale, fill = risposta)) +
  geom_col(width = 0.6, colour = "white") +
  geom_text(aes(label = sprintf("%.1f%%", percentuale)),
            vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = c("1" = "#2166ac", "0" = "#d73027",
                                "NC" = "#a6d96a", "NS" = "#fdae61"),
                    labels = c("1" = "Sì (1)", "0" = "No (0)",
                               "NC" = "Non Class.", "NS" = "Non Sa")) +
  labs(title = "Distribuzione risposte – tutte le categorie",
       subtitle = sprintf("N totale = %d", nrow(dati)),
       x = "Risposta", y = "Percentuale (%)", fill = NULL) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")

print(p1)


# 5b. Grafico proporzioni valide (0/1) per gruppo con IC 95%
ci_gruppi <- dati_binari %>%
  group_by(gruppo) %>%
  summarise(
    n     = n(),
    x     = sum(risposta_num),
    prop  = x / n,
    lower = prop.test(x, n)$conf.int[1],
    upper = prop.test(x, n)$conf.int[2],
    .groups = "drop"
  )

p2 <- ggplot(ci_gruppi, aes(x = gruppo, y = prop, colour = gruppo)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.15, linewidth = 1) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(title = "Proporzione risposte '1' per gruppo (IC 95%)",
       subtitle = "Solo rispondenti validi (esclusi NC e NS)",
       x = "Gruppo", y = "Proporzione", colour = NULL) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")

print(p2)

cat("Script completato.\n")
