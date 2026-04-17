# Summarize 10 columns by response (Yes=1 / No=0) with t-tests
# ---------------------------------------------------------------

set.seed(42)
n <- 200

# --- 1. Simulate dataset ---------------------------------------------------
data <- data.frame(
  response    = sample(c(0, 1), n, replace = TRUE, prob = c(0.45, 0.55)),
  age         = round(rnorm(n, mean = 45, sd = 12)),
  bmi         = round(rnorm(n, mean = 27, sd = 5), 1),
  score       = round(rnorm(n, mean = 50, sd = 10), 1),
  dose_mg     = round(rnorm(n, mean = 200, sd = 50)),
  tto_days    = round(rexp(n, rate = 1/15)),
  systolic_bp = round(rnorm(n, mean = 130, sd = 18)),
  creatinine  = round(rnorm(n, mean = 1.0, sd = 0.3), 2),
  severity    = sample(1:5, n, replace = TRUE),
  hemoglobin  = round(rnorm(n, mean = 13.5, sd = 1.5), 1),
  follow_up   = round(runif(n, min = 1, max = 24))
)

# --- 2. Continuous columns to summarize ------------------------------------
cols <- c("age", "bmi", "score", "dose_mg", "tto_days",
          "systolic_bp", "creatinine", "severity", "hemoglobin", "follow_up")

# --- 3. Summary statistics by group ----------------------------------------
summary_table <- function(df, var, group_var = "response") {
  g0 <- df[[var]][df[[group_var]] == 0]
  g1 <- df[[var]][df[[group_var]] == 1]

  fmt <- function(x) sprintf("%.2f (%.2f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))

  tt <- t.test(g0, g1, var.equal = FALSE)  # Welch t-test

  data.frame(
    Variable    = var,
    No_mean_sd  = fmt(g0),
    Yes_mean_sd = fmt(g1),
    t_statistic = round(tt$statistic, 3),
    df          = round(tt$parameter, 1),
    p_value     = round(tt$p.value, 4),
    Significant = ifelse(tt$p.value < 0.05, "*", ""),
    stringsAsFactors = FALSE
  )
}

results <- do.call(rbind, lapply(cols, summary_table, df = data))
rownames(results) <- NULL

# --- 4. Print results -------------------------------------------------------
cat("=====================================================================\n")
cat("  Summary & t-tests: Yes (1) vs No (0) responses\n")
cat("  Format: Mean (SD)   |   Welch two-sample t-test\n")
cat("=====================================================================\n\n")

cat(sprintf("  %-14s %-20s %-20s %8s %7s %8s\n",
            "Variable", "No (0)", "Yes (1)", "t-stat", "df", "p-value"))
cat(paste(rep("-", 80), collapse = ""), "\n")

for (i in seq_len(nrow(results))) {
  r <- results[i, ]
  cat(sprintf("  %-14s %-20s %-20s %8s %7s %8s %s\n",
              r$Variable, r$No_mean_sd, r$Yes_mean_sd,
              r$t_statistic, r$df, r$p_value, r$Significant))
}

cat(paste(rep("-", 80), collapse = ""), "\n")
cat("  * p < 0.05\n\n")

# --- 5. Group counts -------------------------------------------------------
cat(sprintf("  Group sizes:  No (0) = %d   |   Yes (1) = %d\n\n",
            sum(data$response == 0), sum(data$response == 1)))

# --- 6. Save results to CSV ------------------------------------------------
write.csv(results, "summary_ttest_results.csv", row.names = FALSE)
cat("  Results saved to: summary_ttest_results.csv\n")
