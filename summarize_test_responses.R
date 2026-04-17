# Summarize 10 binary (0/1) columns by response (Yes=1 / No=0) with chi-square
# ------------------------------------------------------------------------------

set.seed(42)
n <- 200

# --- 1. Simulate dataset (all columns are binary 0/1) -------------------------
data <- data.frame(
  response  = sample(c(0, 1), n, replace = TRUE, prob = c(0.45, 0.55)),
  symptom1  = sample(c(0, 1), n, replace = TRUE, prob = c(0.6, 0.4)),
  symptom2  = sample(c(0, 1), n, replace = TRUE, prob = c(0.5, 0.5)),
  symptom3  = sample(c(0, 1), n, replace = TRUE, prob = c(0.7, 0.3)),
  symptom4  = sample(c(0, 1), n, replace = TRUE, prob = c(0.4, 0.6)),
  symptom5  = sample(c(0, 1), n, replace = TRUE, prob = c(0.55, 0.45)),
  symptom6  = sample(c(0, 1), n, replace = TRUE, prob = c(0.65, 0.35)),
  symptom7  = sample(c(0, 1), n, replace = TRUE, prob = c(0.5, 0.5)),
  symptom8  = sample(c(0, 1), n, replace = TRUE, prob = c(0.45, 0.55)),
  symptom9  = sample(c(0, 1), n, replace = TRUE, prob = c(0.6, 0.4)),
  symptom10 = sample(c(0, 1), n, replace = TRUE, prob = c(0.35, 0.65))
)

# --- 2. Columns to summarize --------------------------------------------------
cols <- c("symptom1", "symptom2", "symptom3", "symptom4", "symptom5",
          "symptom6", "symptom7", "symptom8", "symptom9", "symptom10")

# --- 3. Chi-square test per column --------------------------------------------
chisq_table <- function(df, var, group_var = "response") {
  tbl <- table(df[[group_var]], df[[var]])

  n0  <- sum(df[[group_var]] == 0)
  n1  <- sum(df[[group_var]] == 1)
  n0_yes <- tbl["0", "1"]   # column=1 within response=0
  n1_yes <- tbl["1", "1"]   # column=1 within response=1

  # Use Fisher's exact if any expected count < 5
  expected <- chisq.test(tbl)$expected
  test <- if (any(expected < 5)) {
    fisher.test(tbl)
  } else {
    chisq.test(tbl, correct = FALSE)
  }

  method  <- if (inherits(test, "htest") && grepl("Fisher", test$method)) "Fisher" else "Chi-sq"
  stat    <- if (method == "Chi-sq") round(test$statistic, 3) else NA

  data.frame(
    Variable    = var,
    No_n_pct    = sprintf("%d (%.1f%%)", n0_yes, 100 * n0_yes / n0),
    Yes_n_pct   = sprintf("%d (%.1f%%)", n1_yes, 100 * n1_yes / n1),
    Test        = method,
    Statistic   = ifelse(is.na(stat), "-", as.character(stat)),
    p_value     = round(test$p.value, 4),
    Significant = ifelse(test$p.value < 0.05, "*", ""),
    stringsAsFactors = FALSE
  )
}

results <- do.call(rbind, lapply(cols, chisq_table, df = data))
rownames(results) <- NULL

# --- 4. Print results ---------------------------------------------------------
n0 <- sum(data$response == 0)
n1 <- sum(data$response == 1)

cat("=======================================================================\n")
cat("  Summary & Chi-square tests: Yes (1) vs No (0) responses\n")
cat("  Format: n (%)   |   Chi-square or Fisher's exact test\n")
cat("=======================================================================\n\n")
cat(sprintf("  Group sizes:  No (0) = %d   |   Yes (1) = %d   |   Total = %d\n\n",
            n0, n1, n))

cat(sprintf("  %-12s %-16s %-16s %-8s %-9s %-8s\n",
            "Variable", sprintf("No (0) n=%d", n0), sprintf("Yes (1) n=%d", n1),
            "Test", "Statistic", "p-value"))
cat(paste(rep("-", 74), collapse = ""), "\n")

for (i in seq_len(nrow(results))) {
  r <- results[i, ]
  cat(sprintf("  %-12s %-16s %-16s %-8s %-9s %-8s %s\n",
              r$Variable, r$No_n_pct, r$Yes_n_pct,
              r$Test, r$Statistic, r$p_value, r$Significant))
}

cat(paste(rep("-", 74), collapse = ""), "\n")
cat("  * p < 0.05\n\n")

# --- 5. Save to CSV -----------------------------------------------------------
write.csv(results, "summary_chisq_results.csv", row.names = FALSE)
cat("  Results saved to: summary_chisq_results.csv\n")
