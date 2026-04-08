# =============================================================================
# sample_size_rct.R
# -----------------
# Calcolo del campione – Test ONE-SAMPLE a due code sulla proporzione
#
# ENDPOINT PRIMARIO
#   Proporzione di pazienti con >= 5 side effects nei 90 giorni
#   successivi all'ultima somministrazione del farmaco.
#
# DISEGNO
#   Un solo gruppo.  X ~ Binomiale(n, p)
#   H0: p  = 0.30   (soglia di riferimento, pooled storico)
#   H1: p ≠ 0.30   (two-sided)
#
# FORMULA (approssimazione normale alla Binomiale)
#   n = [z_{α/2}·√(p0·(1-p0))  +  z_β·√(p1·(1-p1))]²  /  (p1-p0)²
#
# PARAMETRI: Power=90% | Alpha=5% two-sided | Dropout=10%
# =============================================================================

# install.packages(c("pwr", "ggplot2"))
library(pwr)
library(ggplot2)


# ── 0. Parametri ──────────────────────────────────────────────────────────────

P0      <- 0.30   # soglia H0
ALPHA   <- 0.05   # two-sided
POWER   <- 0.90
DROPOUT <- 0.10

z_a2 <- qnorm(1 - ALPHA / 2)   # 1.9600
z_b  <- qnorm(POWER)            # 1.2816

adj  <- function(n) ceiling(n / (1 - DROPOUT))


# ── 1. Formula analitica (approx. normale alla Binomiale) ─────────────────────

n_formula <- function(p0, p1) {
  d   <- abs(p1 - p0)
  v0  <- p0 * (1 - p0)   # Var Binomiale H0
  v1  <- p1 * (1 - p1)   # Var Binomiale H1
  ceiling(((z_a2 * sqrt(v0) + z_b * sqrt(v1)) / d)^2)
}

# Potenza effettiva (appross. normale) a dato n
power_norm <- function(p0, p1, n) {
  d <- abs(p1 - p0); v1 <- p1 * (1 - p1)
  se <- sqrt(v1 / n)
  pnorm(d / se - z_a2) + pnorm(-d / se - z_a2)
}


# ── 2. Test esatto Binomiale (ricerca n minimo) ───────────────────────────────

# CDF esatta
bcdf <- function(k, n, p) pbinom(k, n, p)

# Potenza esatta two-sided
exact_power <- function(n, p0, p1, alpha) {
  a2   <- alpha / 2
  # soglia alta: min k t.c. P(X >= k | p0) <= alpha/2
  k_hi <- n
  for (k in 1:n) {
    if (1 - bcdf(k - 1, n, p0) <= a2) { k_hi <- k; break }
  }
  # soglia bassa: max k t.c. P(X <= k | p0) <= alpha/2
  k_lo <- -1
  for (k in (n - 1):0) {
    if (bcdf(k, n, p0) <= a2) { k_lo <- k; break }
  }
  pw <- (1 - bcdf(k_hi - 1, n, p1)) + ifelse(k_lo >= 0, bcdf(k_lo, n, p1), 0)
  list(power = pw, k_lo = k_lo, k_hi = k_hi)
}

find_n <- function(p0, p1, alpha, target) {
  for (n in 1:500) {
    res <- exact_power(n, p0, p1, alpha)
    if (res$power >= target)
      return(list(n = n, power = res$power, k_lo = res$k_lo, k_hi = res$k_hi))
  }
}


# ── 3. Calcolo per vari p1 ────────────────────────────────────────────────────

p1_vals <- c(0.40, 0.44, 0.45, 0.50, 0.55, 0.60)

results <- lapply(p1_vals, function(p1) {
  nf   <- n_formula(P0, p1)
  ex   <- find_n(P0, p1, ALPHA, POWER)
  pw_n <- power_norm(P0, p1, nf)
  data.frame(
    p0         = P0,
    p1         = p1,
    delta_pp   = round(abs(p1 - P0) * 100, 0),
    n_formula  = nf,
    n_esatto   = ex$n,
    n_adj      = adj(ex$n),
    power_ex   = round(ex$power * 100, 1),
    k_lo       = ex$k_lo,
    k_hi       = ex$k_hi,
    regola     = paste0("X<=", ex$k_lo, " o X>=", ex$k_hi),
    storico    = p1 == 0.44,
    ufficial   = p1 == 0.50
  )
})

res_df <- do.call(rbind, results)


# ── 4. Stampa risultati ───────────────────────────────────────────────────────

cat("=============================================================\n")
cat(" ONE-SAMPLE BINOMIALE – TEST A DUE CODE\n")
cat(sprintf(" H0: p = %.2f   H1: p ≠ %.2f\n", P0, P0))
cat(sprintf(" Alpha = %.2f two-sided | Power = %.2f | Dropout = %.0f%%\n",
            ALPHA, POWER, DROPOUT * 100))
cat("=============================================================\n\n")

cat(sprintf(" z_{α/2} = %.4f   z_β = %.4f\n\n", z_a2, z_b))

cat(sprintf(" %-6s  %-6s  %-10s  %-9s  %-6s  %-8s  %-22s\n",
            "p1", "Δ(pp)", "n Formula", "n Esatto", "n+DO", "Power", "Regola"))
cat(paste(rep("─", 72), collapse = ""), "\n")

for (i in seq_len(nrow(res_df))) {
  r    <- res_df[i, ]
  flag <- ifelse(r$storico, "  ◄ storico",
          ifelse(r$ufficial, "  ◄ ~n=60 uff.", ""))
  cat(sprintf(" %.2f    %3.0fpp  %10d  %9d  %6d  %7.1f%%  %-22s%s\n",
              r$p1, r$delta_pp, r$n_formula, r$n_esatto,
              r$n_adj, r$power_ex, r$regola, flag))
}

# Dettaglio p1=0.44 e p1=0.50
for (p1_det in c(0.44, 0.50)) {
  r <- res_df[abs(res_df$p1 - p1_det) < 1e-9, ]
  cat("\n──────────────────────────────────────────────────────\n")
  cat(sprintf(" DETTAGLIO p1 = %.2f\n", p1_det))
  cat(sprintf("  n evaluabili    : %d\n", r$n_esatto))
  cat(sprintf("  n + dropout 10%% : %d\n", r$n_adj))
  cat(sprintf("  Power esatta    : %.1f%%\n", r$power_ex))
  cat(sprintf("  Regola          : rifiuta H0 se X <= %d oppure X >= %d\n",
              r$k_lo, r$k_hi))
  cat(sprintf("                    (su %d pazienti evaluabili)\n", r$n_esatto))
}


# ── 5. Verifica con pwr::pwr.p.test (Cohen's h) ──────────────────────────────

cat("\n=============================================================\n")
cat(" VERIFICA – pwr::pwr.p.test (Cohen's h, appross. arcseno)\n")
cat("=============================================================\n")

for (p1 in c(0.44, 0.50)) {
  h  <- ES.h(p1, P0)
  r  <- pwr.p.test(h = abs(h), sig.level = ALPHA, power = POWER,
                   alternative = "two.sided")
  na <- adj(ceiling(r$n))
  cat(sprintf(" p1=%.2f: Cohen h=%.4f  n=%d  n+dropout=%d\n",
              p1, abs(h), ceiling(r$n), na))
}


# ── 6. Curva di potenza ───────────────────────────────────────────────────────

n_seq   <- 20:250
p1_plot <- c(0.44, 0.50, 0.55)

curve_data <- do.call(rbind, lapply(p1_plot, function(p1) {
  data.frame(
    n     = n_seq,
    power = sapply(n_seq, function(n) power_norm(P0, p1, n)),
    p1    = factor(paste0("p1 = ", p1))
  )
}))

gg <- ggplot(curve_data, aes(x = n, y = power * 100,
                              colour = p1, linetype = p1)) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = POWER * 100, linetype = "dashed",
             colour = "black") +
  annotate("text", x = 240, y = POWER * 100 + 1.5,
           label = "Power 90%", hjust = 1, size = 3.2) +
  scale_y_continuous(limits = c(50, 100), breaks = seq(50, 100, 10),
                     labels = function(x) paste0(x, "%")) +
  scale_colour_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A")) +
  labs(
    title    = "Curva di Potenza – ONE-SAMPLE Binomiale",
    subtitle = sprintf("H0: p=%.2f | Alpha=%.0f%% two-sided | X ~ Bin(n, p)",
                       P0, ALPHA * 100),
    x        = "N pazienti (evaluabili)",
    y        = "Potenza (%)",
    colour   = "Alternativa",
    linetype = "Alternativa"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "inside", legend.position.inside = c(0.82, 0.25))

print(gg)


# ── 7. Riepilogo finale ───────────────────────────────────────────────────────

cat("\n=============================================================\n")
cat(" RIEPILOGO FINALE\n")
cat(sprintf(" X~Bin(n,p) | H0: p=%.2f | H1: p≠%.2f | α=%.0f%% | Power=%.0f%% | DO=%.0f%%\n",
            P0, P0, ALPHA*100, POWER*100, DROPOUT*100))
cat("=============================================================\n")
cat(sprintf(" %-8s  %-18s  %-18s\n",
            "p1", "n pazienti enrollati", "Regola decisionale"))
cat(paste(rep("─", 50), collapse=""), "\n")
for (i in seq_len(nrow(res_df))) {
  r <- res_df[i, ]
  cat(sprintf(" %.2f      %6d                   %s\n",
              r$p1, r$n_adj, r$regola))
}
