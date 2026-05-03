# Sommario variabili c2.x con valori: 0, 1, nc, ns
# nc (non compilato) e ns (non so) vengono rimossi
# Confronto: 0 vs 1 con chi-quadro o Fisher esatto
# ------------------------------------------------------------------------------

set.seed(42)
n <- 300

# --- 1. Carica i tuoi dati ----------------------------------------------------
# Sostituisci con:
#   data <- read.csv("tuo_file.csv")
#   data <- readxl::read_excel("tuo_file.xlsx")

vals <- c("0", "1", "nc", "ns")

data <- data.frame(
  response = sample(vals, n, replace = TRUE, prob = c(0.35, 0.50, 0.08, 0.07)),
  c2.2.  = sample(vals, n, replace = TRUE, prob = c(0.30, 0.55, 0.08, 0.07)),
  c2.3.  = sample(vals, n, replace = TRUE, prob = c(0.45, 0.40, 0.08, 0.07)),
  c2.5.  = sample(vals, n, replace = TRUE, prob = c(0.25, 0.60, 0.08, 0.07)),
  c2.7.  = sample(vals, n, replace = TRUE, prob = c(0.50, 0.35, 0.08, 0.07)),
  c2.9.  = sample(vals, n, replace = TRUE, prob = c(0.35, 0.50, 0.08, 0.07)),
  c2.11. = sample(vals, n, replace = TRUE, prob = c(0.40, 0.45, 0.08, 0.07)),
  c2.13. = sample(vals, n, replace = TRUE, prob = c(0.20, 0.65, 0.08, 0.07)),
  c2.15. = sample(vals, n, replace = TRUE, prob = c(0.55, 0.30, 0.08, 0.07)),
  c2.17. = sample(vals, n, replace = TRUE, prob = c(0.30, 0.55, 0.08, 0.07)),
  c2.19. = sample(vals, n, replace = TRUE, prob = c(0.43, 0.42, 0.08, 0.07)),
  c2.21. = sample(vals, n, replace = TRUE, prob = c(0.25, 0.60, 0.08, 0.07)),
  c2.23. = sample(vals, n, replace = TRUE, prob = c(0.47, 0.38, 0.08, 0.07)),
  c2.25. = sample(vals, n, replace = TRUE, prob = c(0.33, 0.52, 0.08, 0.07)),
  c2.27. = sample(vals, n, replace = TRUE, prob = c(0.41, 0.44, 0.08, 0.07)),
  stringsAsFactors = FALSE
)

# --- 2. Colonne da analizzare -------------------------------------------------
cols <- c("c2.2.", "c2.3.", "c2.5.", "c2.7.", "c2.9.", "c2.11.", "c2.13.",
          "c2.15.", "c2.17.", "c2.19.", "c2.21.", "c2.23.", "c2.25.", "c2.27.")

# --- 3. PASSO 1: Frequenze 0 / 1 / nc / ns per ogni variabile ----------------
freq_table <- function(df, var) {
  x     <- as.character(df[[var]])
  total <- length(x)
  data.frame(
    Variabile = var,
    N_1       = sum(x == "1",  na.rm = TRUE),
    Pct_1     = round(100 * sum(x == "1",  na.rm = TRUE) / total, 1),
    N_0       = sum(x == "0",  na.rm = TRUE),
    Pct_0     = round(100 * sum(x == "0",  na.rm = TRUE) / total, 1),
    N_nc      = sum(x == "nc", na.rm = TRUE),
    Pct_nc    = round(100 * sum(x == "nc", na.rm = TRUE) / total, 1),
    N_ns      = sum(x == "ns", na.rm = TRUE),
    Pct_ns    = round(100 * sum(x == "ns", na.rm = TRUE) / total, 1),
    Totale    = total,
    stringsAsFactors = FALSE
  )
}

freq_res <- do.call(rbind, lapply(c("response", cols), freq_table, df = data))
rownames(freq_res) <- NULL

cat("==========================================================================\n")
cat("  PASSO 1 — Frequenze: 1 / 0 / nc / ns per ogni variabile\n")
cat("==========================================================================\n\n")
cat(sprintf("  %-10s  %5s %6s  %5s %6s  %5s %6s  %5s %6s  %7s\n",
            "Variabile", "N 1", "% 1", "N 0", "% 0",
            "N nc", "% nc", "N ns", "% ns", "Totale"))
cat(paste(rep("-", 80), collapse = ""), "\n")
for (i in seq_len(nrow(freq_res))) {
  r <- freq_res[i, ]
  cat(sprintf("  %-10s  %5d %5.1f%%  %5d %5.1f%%  %5d %5.1f%%  %5d %5.1f%%  %7d\n",
              r$Variabile,
              r$N_1, r$Pct_1, r$N_0, r$Pct_0,
              r$N_nc, r$Pct_nc, r$N_ns, r$Pct_ns,
              r$Totale))
}
cat(paste(rep("-", 80), collapse = ""), "\n\n")
write.csv(freq_res, "frequenze_0_1_nc_ns.csv", row.names = FALSE)
cat("  Frequenze salvate in: frequenze_0_1_nc_ns.csv\n\n")

# --- 4. PASSO 2: Chi-quadro 0 vs 1 (esclusi nc e ns) -------------------------
analyze_col <- function(df, var, group_var = "response") {
  # Tieni solo righe con 0 o 1 in entrambe le variabili
  keep <- df[[group_var]] %in% c("0", "1") & df[[var]] %in% c("0", "1")
  sub  <- df[keep, ]
  n_excluded <- sum(!keep)

  g_var   <- as.character(sub[[group_var]])
  col_var <- as.character(sub[[var]])

  n_0 <- sum(g_var == "0")
  n_1 <- sum(g_var == "1")

  # Quanti hanno valore 1 nella colonna per ogni gruppo risposta
  v1_in_g0 <- sum(col_var == "1" & g_var == "0")
  v1_in_g1 <- sum(col_var == "1" & g_var == "1")

  tbl      <- table(g_var, col_var)
  expected <- suppressWarnings(chisq.test(tbl)$expected)
  if (any(expected < 5)) {
    test <- fisher.test(tbl); method <- "Fisher"; stat <- "-"
  } else {
    test <- chisq.test(tbl, correct = FALSE)
    method <- "Chi-sq"; stat <- as.character(round(test$statistic, 3))
  }

  data.frame(
    Variabile     = var,
    N_esclusi     = n_excluded,
    Gr0_1_pct     = sprintf("%d / %d (%.1f%%)", v1_in_g0, n_0, 100 * v1_in_g0 / max(n_0, 1)),
    Gr1_1_pct     = sprintf("%d / %d (%.1f%%)", v1_in_g1, n_1, 100 * v1_in_g1 / max(n_1, 1)),
    Test          = method,
    Statistica    = stat,
    p_value       = round(test$p.value, 4),
    Significativo = ifelse(test$p.value < 0.05, "*", ""),
    stringsAsFactors = FALSE
  )
}

results <- do.call(rbind, lapply(cols, analyze_col, df = data))
rownames(results) <- NULL

data_valid <- data[data$response %in% c("0", "1"), ]
n_0  <- sum(data_valid$response == "0")
n_1  <- sum(data_valid$response == "1")
n_ex <- sum(!data$response %in% c("0", "1"))

cat("==========================================================================\n")
cat("  PASSO 2 — Test: risposta 0 vs 1  (nc e ns esclusi)\n")
cat("  Formato: n con valore=1 / N totale gruppo (%)  |  test per riga\n")
cat("==========================================================================\n\n")
cat(sprintf("  Risposta 0 = %d  |  Risposta 1 = %d  |  nc/ns esclusi = %d\n\n",
            n_0, n_1, n_ex))

cat(sprintf("  %-10s  %-8s  %-22s  %-22s  %-8s  %-9s  %-8s\n",
            "Variabile", "Esclusi",
            sprintf("Risp.0 (n=%d)", n_0),
            sprintf("Risp.1 (n=%d)", n_1),
            "Test", "Statistica", "p-value"))
cat(paste(rep("-", 95), collapse = ""), "\n")
for (i in seq_len(nrow(results))) {
  r <- results[i, ]
  cat(sprintf("  %-10s  %-8s  %-22s  %-22s  %-8s  %-9s  %-8s %s\n",
              r$Variabile, r$N_esclusi,
              r$Gr0_1_pct, r$Gr1_1_pct,
              r$Test, r$Statistica, r$p_value, r$Significativo))
}
cat(paste(rep("-", 95), collapse = ""), "\n")
cat("  * p < 0.05\n\n")

write.csv(results, "summary_chisq_results.csv", row.names = FALSE)
cat("  Risultati salvati in: summary_chisq_results.csv\n")
