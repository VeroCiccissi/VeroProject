# =============================================================================
# sample_size_rct.R
# -----------------
# Calcolo del campione per RCT a due gruppi – Distribuzione Binomiale
#
# Endpoint primario:
#   Proporzione di pazienti con >= 5 side effects (SE) nei 90 giorni
#   successivi all'ultima somministrazione del farmaco.
#
# Modello:
#   X_ctrl  ~ Bin(n, p1)    p1 = 0.44   (dati storici)
#   X_treat ~ Bin(n, p2)    p2 = 0.21   (dati storici)
#   Proporzione ponderata pooled: p = 0.30
#
# Parametri: Power=90%  |  Alpha=5% (two-sided)  |  Dropout=10%
# =============================================================================

# --- Pacchetti necessari -----------------------------------------------------
# install.packages(c("pwr", "ggplot2"))
library(pwr)      # power analysis (Cohen's h, arcsine transform)
library(ggplot2)  # power curve


# =============================================================================
# 0. PARAMETRI GLOBALI
# =============================================================================

P_CTRL   <- 0.44   # X_ctrl  ~ Bin(n, p1)
P_TREAT  <- 0.21   # X_treat ~ Bin(n, p2)
P_POOL   <- 0.30   # proporzione ponderata pooled (media pesata storica)
ALPHA    <- 0.05   # livello di significatività (two-sided)
POWER    <- 0.90   # potenza desiderata
DROPOUT  <- 0.10   # tasso di dropout
MIN_SE   <- 5      # soglia evento: >= 5 SE
FOLLOWUP <- 90     # giorni di follow-up post ultima dose


# =============================================================================
# 1. VERIFICA PROPORZIONE PONDERATA
# =============================================================================

w_ctrl  <- (P_POOL - P_TREAT) / (P_CTRL - P_TREAT)   # peso implicito ctrl
w_treat <- 1 - w_ctrl
p_check <- w_ctrl * P_CTRL + w_treat * P_TREAT

cat("────────────────────────────────────────────────────────\n")
cat("VERIFICA PROPORZIONE PONDERATA (studio storico)\n")
cat("────────────────────────────────────────────────────────\n")
cat(sprintf("  X_ctrl  ~ Bin(n, %.2f)   peso: %.1f%%\n", P_CTRL,  w_ctrl*100))
cat(sprintf("  X_treat ~ Bin(n, %.2f)   peso: %.1f%%\n", P_TREAT, w_treat*100))
cat(sprintf("  p_pool calcolata: %.4f  |  p_pool dichiarata: %.4f\n",
            p_check, P_POOL))
cat("────────────────────────────────────────────────────────\n\n")


# =============================================================================
# 2. ANALISI PRIMARIA – FORMULA FLEISS (varianza Binomiale esplicita)
# =============================================================================
# Formula:
#   n = [z_{α/2}·√(2p̄(1-p̄))  +  z_β·√(p1(1-p1)+p2(1-p2))]²
#       ──────────────────────────────────────────────────────
#                          (p1 − p2)²
#
# La varianza Var[p̂] = p(1-p)/n viene direttamente dalla Binomiale.
# =============================================================================

fleiss_n <- function(p1, p2, alpha = 0.05, power = 0.90) {
  z_a2  <- qnorm(1 - alpha / 2)          # z_{α/2}
  z_b   <- qnorm(power)                  # z_{β}
  p_bar <- (p1 + p2) / 2                 # pooled sotto H0 (alloc. 1:1)
  v1    <- p1 * (1 - p1)                 # Var Binomiale ctrl
  v2    <- p2 * (1 - p2)                 # Var Binomiale treat
  delta <- abs(p1 - p2)
  n_raw <- ((z_a2 * sqrt(2 * p_bar * (1 - p_bar)) +
             z_b  * sqrt(v1 + v2)) / delta)^2
  ceiling(n_raw)
}

# Correzione per continuità (Fleiss 1981) – raccomandata se n < 100
fleiss_cc_n <- function(n, delta) {
  ceiling((n / 4) * (1 + sqrt(1 + 4 / (n * delta)))^2)
}

# Potenza effettiva a n calcolato (approssimazione normale alla Binomiale)
binom_power <- function(p1, p2, n, alpha = 0.05) {
  z_a2  <- qnorm(1 - alpha / 2)
  delta <- abs(p1 - p2)
  se_h1 <- sqrt((p1*(1-p1) + p2*(1-p2)) / n)
  pnorm(delta / se_h1 - z_a2)
}

# --- Calcolo primario --------------------------------------------------------
n_fl    <- fleiss_n(P_CTRL, P_TREAT, ALPHA, POWER)
n_cc    <- fleiss_cc_n(n_fl, abs(P_CTRL - P_TREAT))
n_adj   <- ceiling(n_fl / (1 - DROPOUT))
n_tot   <- n_adj * 2
pw_eff  <- binom_power(P_CTRL, P_TREAT, n_fl, ALPHA)

cat("=============================================================\n")
cat("ANALISI PRIMARIA – N da distribuzione Binomiale (Fleiss)\n")
cat(sprintf("  Endpoint: >= %d SE in %d gg post ultima dose\n", MIN_SE, FOLLOWUP))
cat(sprintf("  X_ctrl~Bin(n,%.2f)  vs  X_treat~Bin(n,%.2f)\n", P_CTRL, P_TREAT))
cat("=============================================================\n")
cat(sprintf("  z_{α/2} = %.4f    z_β = %.4f\n",
            qnorm(1-ALPHA/2), qnorm(POWER)))
cat(sprintf("  Var Bin ctrl  p1(1-p1) = %.4f\n", P_CTRL*(1-P_CTRL)))
cat(sprintf("  Var Bin treat p2(1-p2) = %.4f\n", P_TREAT*(1-P_TREAT)))
cat(sprintf("  Differenza |p1-p2|     = %.2f\n", abs(P_CTRL-P_TREAT)))
cat("-------------------------------------------------------------\n")
cat(sprintf("  N/gruppo  (Fleiss)           : %d\n", n_fl))
cat(sprintf("  N/gruppo  (cont. correction) : %d\n", n_cc))
cat(sprintf("  Potenza effettiva a n_fl     : %.1f%%\n", pw_eff*100))
cat(sprintf("  N/gruppo  + dropout 10%%      : %d\n", n_adj))
cat(sprintf("  N TOTALE  + dropout 10%%      : %d\n", n_tot))
cat("=============================================================\n\n")


# =============================================================================
# 3. VERIFICA – power.prop.test() e pwr::pwr.2p.test()
# =============================================================================

# Base R – usa la stessa approssimazione normale alla Binomiale
res_base <- power.prop.test(
  p1          = P_CTRL,
  p2          = P_TREAT,
  sig.level   = ALPHA,
  power       = POWER,
  alternative = "two.sided"
)

# pwr – usa l'effetto di Cohen h = 2·asin(√p1) − 2·asin(√p2)
h_cohen  <- ES.h(P_CTRL, P_TREAT)          # Cohen's h
res_pwr  <- pwr.2p.test(
  h         = h_cohen,
  sig.level = ALPHA,
  power     = POWER,
  alternative = "two.sided"
)

n_base_r    <- ceiling(res_base$n)
n_base_pwr  <- ceiling(res_pwr$n)
n_adj_base  <- ceiling(n_base_r   / (1 - DROPOUT))
n_adj_pwr   <- ceiling(n_base_pwr / (1 - DROPOUT))

cat("=============================================================\n")
cat("VERIFICA – Funzioni R standard\n")
cat("=============================================================\n")
cat(sprintf("  power.prop.test()  N/gruppo = %d   N totale (+dropout) = %d\n",
            n_base_r, n_adj_base * 2))
cat(sprintf("  pwr.2p.test()      N/gruppo = %d   N totale (+dropout) = %d\n",
            n_base_pwr, n_adj_pwr * 2))
cat(sprintf("  Formula Fleiss     N/gruppo = %d   N totale (+dropout) = %d\n",
            n_fl, n_tot))
cat(sprintf("  [Cohen's h = %.4f]\n", h_cohen))
cat("=============================================================\n\n")


# =============================================================================
# 4. ANALISI DI SENSIBILITÀ – variazione p_treat, p_ctrl fisso a 0.44
# =============================================================================

p_treat_vals <- c(0.15, 0.21, 0.25, 0.30, 0.35)

sens <- lapply(p_treat_vals, function(p2) {
  n_b   <- fleiss_n(P_CTRL, p2, ALPHA, POWER)
  n_c   <- fleiss_cc_n(n_b, abs(P_CTRL - p2))
  pw    <- binom_power(P_CTRL, p2, n_b, ALPHA)
  n_a   <- ceiling(n_b / (1 - DROPOUT))
  data.frame(
    p_ctrl    = P_CTRL,
    p_treat   = p2,
    delta_pp  = round(abs(P_CTRL - p2) * 100, 1),
    n_fleiss  = n_b,
    n_cc      = n_c,
    power_eff = round(pw * 100, 1),
    n_adj     = n_a,
    n_totale  = n_a * 2,
    storico   = ifelse(abs(p2 - P_TREAT) < 1e-9, "(*)", "")
  )
})

sens_df <- do.call(rbind, sens)

cat("=============================================================\n")
cat("SENSIBILITÀ – X_ctrl~Bin(n,0.44) fisso | p_treat variabile\n")
cat("(*) = valore storico principale\n")
cat("=============================================================\n")
print(sens_df, row.names = FALSE)
cat("=============================================================\n\n")


# =============================================================================
# 5. SCENARIO ALTERNATIVO – ONE-SAMPLE  X ~ Bin(n, 0.30)
# =============================================================================
# Proporzione unica ponderata p = 0.30 derivata dalla media pesata dei due
# gruppi storici (44% e 21%).
#
# Formula one-sample (approssimazione normale alla Binomiale):
#   n = [z_{α/2}·√(p0(1-p0))  +  z_β·√(p1(1-p1))]²
#       ────────────────────────────────────────────
#                      (p1 − p0)²
# =============================================================================

onesample_n <- function(p1, p0, alpha = 0.05, power = 0.90) {
  z_a2  <- qnorm(1 - alpha / 2)
  z_b   <- qnorm(power)
  v0    <- p0 * (1 - p0)     # Var Binomiale sotto H0
  v1    <- p1 * (1 - p1)     # Var Binomiale sotto H1
  delta <- abs(p1 - p0)
  n_raw <- ((z_a2 * sqrt(v0) + z_b * sqrt(v1)) / delta)^2
  ceiling(n_raw)
}

p0_vals <- c(0.44, 0.40, 0.21, 0.50)
labels  <- c(
  "p=0.30 < p_ctrl(0.44) → test efficacia",
  "p=0.30 < 0.40 → scenario conservativo",
  "p=0.30 > p_treat(0.21) → test sicurezza",
  "p=0.30 < 0.50 → incidenza sotto soglia"
)

os <- lapply(seq_along(p0_vals), function(i) {
  p0  <- p0_vals[i]
  n_b <- onesample_n(P_POOL, p0, ALPHA, POWER)
  # correzione continuità one-sample: +1/(2·n·delta²)
  n_c <- ceiling(n_b + 1 / (2 * n_b * abs(P_POOL - p0)^2))
  n_a <- ceiling(n_b / (1 - DROPOUT))
  # potenza effettiva
  pw  <- pnorm(abs(P_POOL - p0) / sqrt(P_POOL*(1-P_POOL)/n_b) - qnorm(1-ALPHA/2))
  data.frame(
    p1_H1     = P_POOL,
    p0_H0     = p0,
    delta_pp  = round(abs(P_POOL - p0) * 100, 1),
    n_base    = n_b,
    n_cc      = n_c,
    power_eff = round(pw * 100, 1),
    n_adj     = n_a,
    interp    = labels[i]
  )
})

os_df <- do.call(rbind, os)

cat("=============================================================\n")
cat("SCENARIO B – One-sample  X ~ Bin(n, 0.30)\n")
cat(sprintf("  Proporzione ponderata p=%.2f | Power=%.0f%% | Alpha=%.0f%% | Dropout=%.0f%%\n",
            P_POOL, POWER*100, ALPHA*100, DROPOUT*100))
cat("=============================================================\n")
print(os_df, row.names = FALSE)
cat("=============================================================\n\n")

# Verifica con prop.test power base R (one-sample)
cat("VERIFICA one-sample con pwr::pwr.p.test()\n")
cat("─────────────────────────────────────────\n")
for (i in seq_along(p0_vals)) {
  h_os  <- ES.h(P_POOL, p0_vals[i])     # Cohen's h one-sample
  res_os <- pwr.p.test(
    h         = abs(h_os),
    sig.level = ALPHA,
    power     = POWER,
    alternative = "two.sided"
  )
  n_os_adj <- ceiling(res_os$n / (1 - DROPOUT))
  cat(sprintf("  H0: p=%.2f   n=%d   n+dropout=%d\n",
              p0_vals[i], ceiling(res_os$n), n_os_adj))
}
cat("\n")


# =============================================================================
# 6. CONFRONTO DISEGNI
# =============================================================================

confronto <- data.frame(
  disegno = c(
    "Two-sample: X_ctrl~Bin(n,0.44) vs X_treat~Bin(n,0.21)",
    "One-sample: X~Bin(n,0.30)  vs  H0: p=0.44",
    "One-sample: X~Bin(n,0.30)  vs  H0: p=0.40"
  ),
  n_per_gruppo = c(n_fl, NA, NA),
  n_totale_raw = c(n_fl * 2,
                   os_df$n_base[os_df$p0_H0 == 0.44],
                   os_df$n_base[os_df$p0_H0 == 0.40]),
  n_totale_adj = c(n_tot,
                   os_df$n_adj[os_df$p0_H0 == 0.44],
                   os_df$n_adj[os_df$p0_H0 == 0.40])
)

cat("=============================================================\n")
cat("CONFRONTO DISEGNI (stesso Power 90%, Alpha 5%, Dropout 10%)\n")
cat("=============================================================\n")
print(confronto, row.names = FALSE)
cat("=============================================================\n\n")


# =============================================================================
# 7. CURVA DI POTENZA – ggplot2
# =============================================================================

# Griglia n x p_treat per la curva two-sample
n_seq      <- 30:250
p2_vals    <- c(0.35, 0.30, 0.25, 0.21, 0.15)

curve_data <- expand.grid(n = n_seq, p_treat = p2_vals) |>
  transform(
    power = mapply(binom_power, P_CTRL, p_treat, n, ALPHA),
    label = paste0("p_treat=", p_treat)
  )

gg_two <- ggplot(curve_data, aes(x = n, y = power * 100,
                                  colour = label, linetype = label)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = POWER * 100, linetype = "dashed",
             colour = "black", linewidth = 0.6) +
  annotate("text", x = max(n_seq) * 0.98, y = POWER * 100 + 1.5,
           label = "Power 90%", hjust = 1, size = 3.2) +
  scale_y_continuous(limits = c(50, 100), breaks = seq(50, 100, 10),
                     labels = function(x) paste0(x, "%")) +
  scale_colour_brewer(palette = "Set1") +
  labs(
    title    = "Curva di Potenza – RCT Two-sample",
    subtitle = sprintf("X_ctrl~Bin(n, %.2f) fisso | Alpha=%.0f%% | Two-sided",
                       P_CTRL, ALPHA * 100),
    x        = "N per gruppo",
    y        = "Potenza (%)",
    colour   = "P trattamento",
    linetype = "P trattamento"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "inside", legend.position.inside = c(0.82, 0.25))

print(gg_two)


# Curva one-sample (p = 0.30 vs vari p0)
n_seq_os <- 20:300
p0_curve <- c(0.44, 0.40, 0.21, 0.50)

curve_os <- expand.grid(n = n_seq_os, p0 = p0_curve) |>
  transform(
    power = mapply(function(n, p0)
      pnorm(abs(P_POOL - p0) / sqrt(P_POOL*(1-P_POOL)/n) - qnorm(1-ALPHA/2)),
      n, p0),
    label = paste0("H0: p=", p0)
  )

gg_one <- ggplot(curve_os, aes(x = n, y = power * 100,
                                colour = label, linetype = label)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = POWER * 100, linetype = "dashed",
             colour = "black", linewidth = 0.6) +
  annotate("text", x = max(n_seq_os) * 0.98, y = POWER * 100 + 1.5,
           label = "Power 90%", hjust = 1, size = 3.2) +
  scale_y_continuous(limits = c(50, 100), breaks = seq(50, 100, 10),
                     labels = function(x) paste0(x, "%")) +
  scale_colour_brewer(palette = "Dark2") +
  labs(
    title    = "Curva di Potenza – One-sample",
    subtitle = sprintf("X~Bin(n, %.2f) | Alpha=%.0f%% | Two-sided",
                       P_POOL, ALPHA * 100),
    x        = "N totale",
    y        = "Potenza (%)",
    colour   = "Ipotesi nulla",
    linetype = "Ipotesi nulla"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "inside", legend.position.inside = c(0.82, 0.25))

print(gg_one)


# =============================================================================
# 8. RIEPILOGO FINALE
# =============================================================================

cat("\n")
cat("══════════════════════════════════════════════════════════════\n")
cat(" RIEPILOGO – SAMPLE SIZE  (Power=90% | Alpha=5% | Dropout=10%)\n")
cat("══════════════════════════════════════════════════════════════\n")
cat(sprintf(" Endpoint: >= %d SE nei %d gg post ultima dose\n", MIN_SE, FOLLOWUP))
cat("\n")
cat(" [A] TWO-SAMPLE  X_ctrl~Bin(n,0.44) vs X_treat~Bin(n,0.21)\n")
cat(sprintf("     N/gruppo (Fleiss)     : %d\n", n_fl))
cat(sprintf("     N/gruppo (cont.corr.) : %d\n", n_cc))
cat(sprintf("     N TOTALE + dropout    : %d\n", n_tot))
cat("\n")
cat(" [B] ONE-SAMPLE  X~Bin(n,0.30)  [media ponderata]\n")
for (i in seq_along(p0_vals)) {
  cat(sprintf("     H0: p=%.2f   N totale + dropout = %d\n",
              p0_vals[i], os_df$n_adj[i]))
}
cat("══════════════════════════════════════════════════════════════\n")
