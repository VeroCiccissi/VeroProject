# =============================================================================
# sample_size_rct.R
# -----------------
# Calcolo del campione – Test ONE-SAMPLE ONE-SIDED sulla proporzione
#
# ENDPOINT PRIMARIO
#   Proporzione di pazienti con >= 5 side effects nei 90 giorni
#   successivi all'ultima somministrazione del farmaco.
#
# DISEGNO
#   Un solo gruppo.  X ~ Binomiale(n, p)
#   H0: p  = 0.30   (soglia di riferimento, pooled storico)
#   H1: p  > 0.30   (ONE-SIDED)
#
# FORMULA (approssimazione normale, one-sided)
#   n = [z_α · √(p0·(1-p0))  +  z_β · √(p1·(1-p1))]²  /  (p1-p0)²
#   z_α = qnorm(0.95) = 1.6449  (one-sided 5%)
#
# NOTE
#   p1=0.50, one-sided → n_esatto=53 → +dropout=59 ≈ 60 (documento ufficiale)
#   p1=0.44, one-sided → n_esatto=105 → +dropout=117 (scenario storico)
#
# PARAMETRI: Power=90% | Alpha=5% one-sided | Dropout=10%
# =============================================================================

# install.packages(c("pwr", "ggplot2"))
library(pwr)
library(ggplot2)


# ── 0. Parametri ──────────────────────────────────────────────────────────────

P0      <- 0.30
ALPHA   <- 0.05   # one-sided
POWER   <- 0.90
DROPOUT <- 0.10

z_a <- qnorm(1 - ALPHA)    # 1.6449  (one-sided)
z_b <- qnorm(POWER)         # 1.2816

adj <- function(n) ceiling(n / (1 - DROPOUT))


# ── 1. Formula analitica (one-sided) ─────────────────────────────────────────

n_formula <- function(p0, p1) {
  d <- p1 - p0
  ceiling(((z_a * sqrt(p0*(1-p0)) + z_b * sqrt(p1*(1-p1))) / d)^2)
}

power_norm <- function(p0, p1, n) {
  pnorm((p1 - p0) / sqrt(p1*(1-p1)/n) - z_a)
}


# ── 2. Test esatto Binomiale (one-sided) ──────────────────────────────────────

bcdf <- function(k, n, p) pbinom(k, n, p)

exact_power_1s <- function(n, p0, p1, alpha) {
  k_hi <- n
  for (k in 1:n) {
    if (1 - bcdf(k - 1, n, p0) <= alpha) { k_hi <- k; break }
  }
  pw <- 1 - bcdf(k_hi - 1, n, p1)
  list(power = pw, k_hi = k_hi)
}

find_n <- function(p0, p1, alpha, target) {
  for (n in 1:500) {
    res <- exact_power_1s(n, p0, p1, alpha)
    if (res$power >= target)
      return(list(n = n, power = res$power, k_hi = res$k_hi))
  }
}


# ── 3. Calcolo per vari p1 ────────────────────────────────────────────────────

p1_vals <- c(0.40, 0.44, 0.45, 0.50, 0.55, 0.60)

results <- lapply(p1_vals, function(p1) {
  nf  <- n_formula(P0, p1)
  ex  <- find_n(P0, p1, ALPHA, POWER)
  data.frame(
    p0        = P0,
    p1        = p1,
    delta_pp  = round((p1 - P0) * 100),
    n_formula = nf,
    n_esatto  = ex$n,
    n_adj     = adj(ex$n),
    power_ex  = round(ex$power * 100, 1),
    k_hi      = ex$k_hi,
    storico   = p1 == 0.44,
    ufficial  = p1 == 0.50
  )
})

res_df <- do.call(rbind, results)


# ── 4. Output ─────────────────────────────────────────────────────────────────

cat("=============================================================\n")
cat(" ONE-SAMPLE BINOMIALE – TEST ONE-SIDED\n")
cat(sprintf(" H0: p = %.2f   H1: p > %.2f\n", P0, P0))
cat(sprintf(" Alpha = %.2f one-sided | Power = %.2f | Dropout = %.0f%%\n",
            ALPHA, POWER, DROPOUT * 100))
cat(sprintf(" z_alpha(one-sided) = %.4f   z_beta = %.4f\n", z_a, z_b))
cat("=============================================================\n\n")

cat(sprintf(" %-5s  %-6s  %-10s  %-9s  %-6s  %-8s  %-10s\n",
            "p1", "Δ(pp)", "n Formula", "n Esatto", "n+DO", "Power", "Regola"))
cat(paste(rep("─", 62), collapse = ""), "\n")

for (i in seq_len(nrow(res_df))) {
  r    <- res_df[i, ]
  flag <- if (r$storico)  "  ◄ storico (117 pz)" else
          if (r$ufficial) "  ◄ MATCH doc.uff. (~60 pz)" else ""
  cat(sprintf(" %.2f    %3dpp  %10d  %9d  %6d  %7.1f%%  X>=%d%s\n",
              r$p1, r$delta_pp, r$n_formula, r$n_esatto,
              r$n_adj, r$power_ex, r$k_hi, flag))
}


# ── 5. Dettaglio scenari chiave ───────────────────────────────────────────────

cat("\n")
for (p1_det in c(0.44, 0.50)) {
  r <- res_df[abs(res_df$p1 - p1_det) < 1e-9, ]
  cat("──────────────────────────────────────────────────────\n")
  cat(sprintf(" DETTAGLIO  p1 = %.2f\n", p1_det))
  cat(sprintf("  n evaluabili       : %d\n",   r$n_esatto))
  cat(sprintf("  n enrollati (+DO)  : %d\n",   r$n_adj))
  cat(sprintf("  Power esatta       : %.1f%%\n", r$power_ex))
  cat(sprintf("  Regola decisionale : rifiuta H0 se X >= %d\n", r$k_hi))
  cat(sprintf("  Interpretazione    : >= %d/%d = %.1f%% pz con >=5 SE in 90gg\n",
              r$k_hi, r$n_esatto, r$k_hi / r$n_esatto * 100))
}


# ── 6. Verifica con pwr::pwr.p.test (Cohen's h) ──────────────────────────────

cat("\n──────────────────────────────────────────────────────\n")
cat(" VERIFICA pwr::pwr.p.test (Cohen's h, one-sided)\n")
cat("──────────────────────────────────────────────────────\n")
for (p1 in c(0.44, 0.50)) {
  h   <- ES.h(p1, P0)
  res <- pwr.p.test(h = h, sig.level = ALPHA, power = POWER,
                    alternative = "greater")
  cat(sprintf(" p1=%.2f: h=%.4f  n=%d  n+DO=%d\n",
              p1, h, ceiling(res$n), adj(ceiling(res$n))))
}


# ── 7. Curva di potenza ───────────────────────────────────────────────────────

n_seq   <- 20:200
p1_plot <- c(0.44, 0.50, 0.55)

curve_data <- do.call(rbind, lapply(p1_plot, function(p1) {
  data.frame(
    n     = n_seq,
    power = sapply(n_seq, function(n) power_norm(P0, p1, n)),
    grp   = factor(sprintf("p1 = %.2f", p1))
  )
}))

gg <- ggplot(curve_data, aes(x = n, y = power * 100, colour = grp, linetype = grp)) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = POWER * 100, linetype = "dashed", colour = "black") +
  annotate("text", x = 190, y = POWER * 100 + 1.5,
           label = "Power 90%", hjust = 1, size = 3.2) +
  geom_vline(xintercept = 53, linetype = "dotted", colour = "grey40") +
  annotate("text", x = 55, y = 55, label = "n=53\n(p1=0.50)", hjust = 0, size = 3) +
  geom_vline(xintercept = 105, linetype = "dotted", colour = "grey40") +
  annotate("text", x = 107, y = 55, label = "n=105\n(p1=0.44)", hjust = 0, size = 3) +
  scale_y_continuous(limits = c(50, 100), breaks = seq(50, 100, 10),
                     labels = function(x) paste0(x, "%")) +
  scale_colour_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A")) +
  labs(
    title    = "Curva di Potenza – ONE-SAMPLE ONE-SIDED Binomiale",
    subtitle = sprintf("H0: p=%.2f | H1: p>%.2f | Alpha=%.0f%% one-sided",
                       P0, P0, ALPHA * 100),
    x        = "N pazienti enrollati (post dropout)",
    y        = "Potenza (%)",
    colour   = "Alternativa p1",
    linetype = "Alternativa p1"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "inside", legend.position.inside = c(0.82, 0.25))

print(gg)


# ── 8. Riepilogo finale ───────────────────────────────────────────────────────

cat("\n=============================================================\n")
cat(" RIEPILOGO FINALE\n")
cat(sprintf(" X~Bin(n,p) | H0:p=%.2f | H1:p>%.2f | α=%.0f%% one-sided | Power=%.0f%% | DO=%.0f%%\n",
            P0, P0, ALPHA*100, POWER*100, DROPOUT*100))
cat("=============================================================\n")
for (i in seq_len(nrow(res_df))) {
  r    <- res_df[i, ]
  note <- if (r$storico)  " ← scenario storico" else
          if (r$ufficial) " ← documento ufficiale (~60)" else ""
  cat(sprintf(" p1=%.2f:  n=%d pazienti enrollati  |  rifiuta H0 se X>=%d%s\n",
              r$p1, r$n_adj, r$k_hi, note))
}
