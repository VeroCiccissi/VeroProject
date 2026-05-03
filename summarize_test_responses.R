# Sommario variabili c2.x (si/no/non so) confrontate tra si vs no
# I "non so" vengono rimossi prima dell'analisi
# Test: Chi-quadro o Fisher esatto (automatico)
# ------------------------------------------------------------------------------

set.seed(42)
n <- 300

# --- 1. Carica i tuoi dati ----------------------------------------------------
# Sostituisci questa sezione con:
#   data <- read.csv("tuo_file.csv")
#   oppure
#   data <- readxl::read_excel("tuo_file.xlsx")

vals <- c("si", "no", "non so")

data <- data.frame(
  response = sample(vals, n, replace = TRUE, prob = c(0.5, 0.35, 0.15)),
  c2.2.  = sample(vals, n, replace = TRUE, prob = c(0.55, 0.30, 0.15)),
  c2.3.  = sample(vals, n, replace = TRUE, prob = c(0.40, 0.45, 0.15)),
  c2.5.  = sample(vals, n, replace = TRUE, prob = c(0.60, 0.25, 0.15)),
  c2.7.  = sample(vals, n, replace = TRUE, prob = c(0.35, 0.50, 0.15)),
  c2.9.  = sample(vals, n, replace = TRUE, prob = c(0.50, 0.35, 0.15)),
  c2.11. = sample(vals, n, replace = TRUE, prob = c(0.45, 0.40, 0.15)),
  c2.13. = sample(vals, n, replace = TRUE, prob = c(0.65, 0.20, 0.15)),
  c2.15. = sample(vals, n, replace = TRUE, prob = c(0.30, 0.55, 0.15)),
  c2.17. = sample(vals, n, replace = TRUE, prob = c(0.55, 0.30, 0.15)),
  c2.19. = sample(vals, n, replace = TRUE, prob = c(0.42, 0.43, 0.15)),
  c2.21. = sample(vals, n, replace = TRUE, prob = c(0.60, 0.25, 0.15)),
  c2.23. = sample(vals, n, replace = TRUE, prob = c(0.38, 0.47, 0.15)),
  c2.25. = sample(vals, n, replace = TRUE, prob = c(0.52, 0.33, 0.15)),
  c2.27. = sample(vals, n, replace = TRUE, prob = c(0.44, 0.41, 0.15)),
  stringsAsFactors = FALSE
)

# --- 2. Nomi colonne da analizzare --------------------------------------------
cols <- c("c2.2.", "c2.3.", "c2.5.", "c2.7.", "c2.9.", "c2.11.", "c2.13.",
          "c2.15.", "c2.17.", "c2.19.", "c2.21.", "c2.23.", "c2.25.", "c2.27.")

# --- 3. FREQUENZE ASSOLUTE E RELATIVE per ogni variabile (si / no / non so) --
freq_table <- function(df, var) {
  x     <- df[[var]]
  total <- length(x)
  n_si  <- sum(x == "si",     na.rm = TRUE)
  n_no  <- sum(x == "no",     na.rm = TRUE)
  n_ns  <- sum(x == "non so", na.rm = TRUE)
  data.frame(
    Variabile = var,
    N_si      = n_si,
    Pct_si    = round(100 * n_si / total, 1),
    N_no      = n_no,
    Pct_no    = round(100 * n_no / total, 1),
    N_nonso   = n_ns,
    Pct_nonso = round(100 * n_ns / total, 1),
    Totale    = total,
    stringsAsFactors = FALSE
  )
}

# Includi anche la variabile risposta nel conteggio
all_vars   <- c("response", cols)
freq_res   <- do.call(rbind, lapply(all_vars, freq_table, df = data))
rownames(freq_res) <- NULL

cat("==========================================================================\n")
cat("  PASSO 1 — Frequenze: SI / NO / NON SO per ogni variabile\n")
cat("==========================================================================\n\n")
cat(sprintf("  %-10s  %6s %7s  %6s %7s  %8s %9s  %7s\n",
            "Variabile", "N si", "% si", "N no", "% no",
            "N non so", "% non so", "Totale"))
cat(paste(rep("-", 72), collapse = ""), "\n")
for (i in seq_len(nrow(freq_res))) {
  r <- freq_res[i, ]
  cat(sprintf("  %-10s  %6d %6.1f%%  %6d %6.1f%%  %8d %8.1f%%  %7d\n",
              r$Variabile, r$N_si, r$Pct_si,
              r$N_no, r$Pct_no,
              r$N_nonso, r$Pct_nonso, r$Totale))
}
cat(paste(rep("-", 72), collapse = ""), "\n\n")

# Salva frequenze
write.csv(freq_res, "frequenze_si_no_nonso.csv", row.names = FALSE)
cat("  Frequenze salvate in: frequenze_si_no_nonso.csv\n\n")

# --- 4. Funzione analisi per ogni colonna (si vs no, senza non so) ------------
analyze_col <- function(df, var, group_var = "response") {

  # Rimuovi "non so" sia dalla variabile che dalla risposta
  keep <- df[[group_var]] != "non so" & df[[var]] != "non so"
  sub  <- df[keep, ]

  n_removed <- sum(!keep)

  g_var   <- sub[[group_var]]
  col_var <- sub[[var]]

  n_no  <- sum(g_var == "no")
  n_si  <- sum(g_var == "si")

  # n e % di "si" nella colonna per ciascun gruppo risposta
  si_in_no <- sum(col_var == "si" & g_var == "no")
  si_in_si <- sum(col_var == "si" & g_var == "si")

  tbl <- table(g_var, col_var)  # righe = gruppi risposta, colonne = si/no

  # Scegli test: Fisher se atteso < 5, altrimenti chi-quadro
  expected <- suppressWarnings(chisq.test(tbl)$expected)
  if (any(expected < 5)) {
    test   <- fisher.test(tbl)
    method <- "Fisher"
    stat   <- "-"
  } else {
    test   <- chisq.test(tbl, correct = FALSE)
    method <- "Chi-sq"
    stat   <- as.character(round(test$statistic, 3))
  }

  data.frame(
    Variabile      = var,
    N_rimossi      = n_removed,
    No_n_pct       = sprintf("%d / %d (%.1f%%)", si_in_no, n_no, 100 * si_in_no / max(n_no, 1)),
    Si_n_pct       = sprintf("%d / %d (%.1f%%)", si_in_si, n_si, 100 * si_in_si / max(n_si, 1)),
    Test           = method,
    Statistica     = stat,
    p_value        = round(test$p.value, 4),
    Significativo  = ifelse(test$p.value < 0.05, "*", ""),
    stringsAsFactors = FALSE
  )
}

results <- do.call(rbind, lapply(cols, analyze_col, df = data))
rownames(results) <- NULL

# --- 5. Stampa risultati confronto SI vs NO -----------------------------------
# Conta dopo rimozione "non so" dalla variabile risposta
data_clean <- data[data$response != "non so", ]
n_no  <- sum(data_clean$response == "no")
n_si  <- sum(data_clean$response == "si")
n_ns  <- sum(data$response == "non so")

cat("==========================================================================\n")
cat("  PASSO 2 — Test: risposta SI vs NO  (\"non so\" esclusi)\n")
cat("  Formato colonne: n SI / N totale gruppo (%)  |  test per riga\n")
cat("==========================================================================\n\n")
cat(sprintf("  Risposta NO  = %d  |  Risposta SI = %d  |  Non so rimossi = %d\n\n",
            n_no, n_si, n_ns))

cat(sprintf("  %-10s  %-7s  %-22s  %-22s  %-8s  %-9s  %-8s\n",
            "Variabile", "Rimossi",
            sprintf("NO (n=%d)  si/tot(%%)", n_no),
            sprintf("SI (n=%d)  si/tot(%%)", n_si),
            "Test", "Statistica", "p-value"))
cat(paste(rep("-", 95), collapse = ""), "\n")

for (i in seq_len(nrow(results))) {
  r <- results[i, ]
  cat(sprintf("  %-10s  %-7s  %-22s  %-22s  %-8s  %-9s  %-8s %s\n",
              r$Variabile, r$N_rimossi,
              r$No_n_pct, r$Si_n_pct,
              r$Test, r$Statistica, r$p_value, r$Significativo))
}

cat(paste(rep("-", 95), collapse = ""), "\n")
cat("  * p < 0.05\n\n")

# --- 6. Salva CSV confronto ---------------------------------------------------
write.csv(results, "summary_chisq_results.csv", row.names = FALSE)
cat("  Risultati salvati in: summary_chisq_results.csv\n")
