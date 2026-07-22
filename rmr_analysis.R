############################################################
# LIBRERIE
############################################################

library(dplyr)
library(tidyr)
library(stringr)
library(rlang)
library(haven)
library(purrr)

############################################################
# LOAD DATASET
############################################################

ae_cases     <- read_sas("//vivaldi/ShBiostatistics/PV/Argus estrazioni/Validation/Argus2SAS_Passthrough/ae_cases.sas7bdat")
ae_dosages   <- read_sas("//vivaldi/ShBiostatistics/PV/Argus estrazioni/Validation/Argus2SAS_Passthrough/ae_dosages.sas7bdat")
ae_drugs     <- read_sas("//vivaldi/ShBiostatistics/PV/Argus estrazioni/Validation/Argus2SAS_Passthrough/ae_drugs.sas7bdat")
ae_events    <- read_sas("//vivaldi/ShBiostatistics/PV/Argus estrazioni/Validation/Argus2SAS_Passthrough/ae_events.sas7bdat")
date_details <- read_sas("//vivaldi/ShBiostatistics/PV/Argus estrazioni/Validation/Argus2SAS_Passthrough/date_details.sas7bdat")

############################################################
# PARAMETRI TEMPORALI
############################################################

datefrom     <- as.Date("2026-04-01")
dateto       <- as.Date("2026-06-30")
lastdatefrom <- as.Date("2026-01-01")
lastdateto   <- as.Date("2026-03-31")

############################################################
# DATASET BASE CUMULATIVO
############################################################

cumulative <- ae_events %>%
  inner_join(ae_drugs %>% filter(drug_role == "Suspect"), by = "CASE_ID") %>%
  inner_join(ae_cases,     by = "CASE_ID") %>%
  inner_join(date_details, by = c("CASE_ID", "CASE_NUM")) %>%
  filter(INIT_REPT_DATE <= dateto, preferred_event_term != "") %>%
  transmute(
    CASE_ID, CASE_NUM,
    event_seq_num, event_llt, preferred_event_term, body_system, event_did,
    prod_seq_num,
    drug_name    = toupper(drug_name),
    GENERIC_NAME,
    Serious_flag_event, died_flag_event,
    Serious_flag_case,  died_flag_case,
    recdate = INIT_REPT_DATE
  ) %>%
  left_join(
    ae_dosages %>% select(CASE_ID, prod_seq_num, eff_flag),
    by = c("CASE_ID", "prod_seq_num")
  )

############################################################
# DATASET TRIMESTRALI
# (cumulative già filtrato su <= dateto, quindi
#  filtrare >= datefrom dà l'intervallo corretto)
############################################################

interval     <- cumulative %>% filter(recdate >= datefrom)
lastinterval <- cumulative %>% filter(recdate >= lastdatefrom, recdate <= lastdateto)

############################################################
# HELPER – filtro standard "prodotto aziendale"
# Riduce la ripetizione nei 15+ prodotti della lista
############################################################

drug_f <- function(generic, extra = NULL) {
  base <- sprintf(
    "GENERIC_NAME == '%s' & drug_name != '' & !str_detect(drug_name, 'NOT COMPANY')",
    generic
  )
  if (!is.null(extra)) paste0("(", base, ") & (", extra, ")") else base
}

############################################################
# FUNZIONE RMR
############################################################

RMR <- function(drug_filter, active_substance) {

  drug_expr         <- parse_expr(drug_filter)
  cumulative_drug   <- cumulative    %>% filter(!!drug_expr)
  interval_drug     <- interval      %>% filter(!!drug_expr)
  lastinterval_drug <- lastinterval  %>% filter(!!drug_expr)

  n_cum  <- n_distinct(cumulative_drug$CASE_ID)
  n_intv <- n_distinct(interval_drug$CASE_ID)

  # Fix: use || (scalar) instead of | (vectorised) for early-return check
  if (n_cum == 0 || n_intv == 0) {
    return(tibble(ACTIVE_SUBSTANCE = active_substance, RESULT = "NO DATA FOUND"))
  }

  # ── Conteggi base della tabella 2x2 ──────────────────────────────────────
  N   <- n_distinct(cumulative$CASE_ID)        # totale casi nel database
  ApB <- n_cum                                  # totale casi per il farmaco (a+b)

  a <- cumulative_drug %>%
    group_by(body_system, preferred_event_term) %>%
    summarise(A = n_distinct(CASE_ID), .groups = "drop")

  apc <- cumulative %>%
    group_by(body_system, preferred_event_term) %>%
    summarise(ApC = n_distinct(CASE_ID), .groups = "drop")

  fatal_PT <- cumulative_drug %>%
    filter(died_flag_event == "Y" | died_flag_case == "Y") %>%
    group_by(body_system, preferred_event_term) %>%
    summarise(TotalFatal = n_distinct(CASE_ID), .groups = "drop")

  intvcount_PT <- interval_drug %>%
    group_by(body_system, preferred_event_term) %>%
    summarise(New = n_distinct(CASE_ID), .groups = "drop")

  intvfatal_PT <- interval_drug %>%
    filter(died_flag_event == "Y" | died_flag_case == "Y") %>%
    group_by(body_system, preferred_event_term) %>%
    summarise(NewFatal = n_distinct(CASE_ID), .groups = "drop")

  lastcount_PT <- lastinterval_drug %>%
    group_by(body_system, preferred_event_term) %>%
    summarise(Last = n_distinct(CASE_ID), .groups = "drop")

  # ── Metriche (tutte vettorizzate – rimosso rowwise()) ────────────────────
  metrics <- a %>%
    left_join(apc, by = c("body_system", "preferred_event_term")) %>%
    mutate(
      B = ApB - A,
      C = ApC - A,
      D = (N - ApB) - C
    ) %>%
    mutate(
      Observed = A,
      Expected = ((A + B) * (A + C)) / (A + B + C + D),
      OE_Ratio = Observed / Expected,

      # PRR
      PRR     = if_else(C > 0, (A / (A + B)) / (C / (C + D)),   NA_real_),
      PRR_SE  = if_else(C > 0, sqrt(1/A - 1/(A+B) + 1/C - 1/(C+D)), NA_real_),
      PRR_LCL = exp(log(PRR) - 1.96 * PRR_SE),
      PRR_UCL = exp(log(PRR) + 1.96 * PRR_SE),

      # ROR
      ROR     = if_else(B > 0 & C > 0 & D > 0, (A * D) / (B * C),            NA_real_),
      ROR_SE  = if_else(B > 0 & C > 0 & D > 0, sqrt(1/A + 1/B + 1/C + 1/D), NA_real_),
      ROR_LCL = exp(log(ROR) - 1.96 * ROR_SE),
      ROR_UCL = exp(log(ROR) + 1.96 * ROR_SE),

      # Chi-quadro (non corretto per continuità)
      CHISQ = ((A+B+C+D) * (A*D - B*C)^2) / ((A+B) * (C+D) * (A+C) * (B+D)),

      # IC (Information Component, approssimazione Bayesiana semplificata)
      IC    = log2((Observed + 0.5) / (Expected + 0.5)),
      IC_SE = sqrt(1 / (Observed + 0.5)),
      IC025 = IC - 1.96 * IC_SE,
      IC975 = IC + 1.96 * IC_SE,

      # EBGM (rapporto O/E con shrinkage)
      EBGM  = (Observed + 0.5) / (Expected + 0.5)
    ) %>%
    mutate(
      # Fisher exact: unica operazione che richiede iterazione per riga
      FISHER_P = pmap_dbl(
        list(A, B, C, D),
        function(a, b, c, d) {
          tryCatch(
            fisher.test(matrix(c(a, b, c, d), nrow = 2))$p.value,
            error = function(e) NA_real_
          )
        }
      ),

      # Fix: aggiunto !is.na(PRR) – senza il check, MHRA_SIGNAL era NA quando PRR è NA
      MHRA_SIGNAL = if_else(!is.na(PRR) & A >= 3 & PRR >= 2 & CHISQ >= 4, "YES", "NO"),
      ROR_SIGNAL  = if_else(!is.na(ROR_LCL) & ROR_LCL > 1,               "YES", "NO"),
      IC_SIGNAL   = if_else(!is.na(IC025)   & IC025   > 0,               "YES", "NO"),
      EBGM_SIGNAL = if_else(!is.na(EBGM)    & EBGM   >= 2,               "YES", "NO")
    ) %>%
    rename(Total = A)

  # ── Assembla risultati finali ─────────────────────────────────────────────
  # Fix: left_join (non full_join) – le tabelle di destra sono sottoinsiemi
  results <- metrics %>%
    left_join(intvcount_PT,  by = c("body_system", "preferred_event_term")) %>%
    left_join(lastcount_PT,  by = c("body_system", "preferred_event_term")) %>%
    left_join(fatal_PT,      by = c("body_system", "preferred_event_term")) %>%
    left_join(intvfatal_PT,  by = c("body_system", "preferred_event_term")) %>%
    replace_na(list(New = 0L, Last = 0L, TotalFatal = 0L, NewFatal = 0L)) %>%
    filter(New > 0) %>%
    mutate(
      GrowthRatio    = round(New / (Last + 1), 2),
      QoQ_Percent    = round(100 * (New - Last) / pmax(Last, 1), 2),
      FatalityRatio  = round(100 * TotalFatal / Total, 2),
      SignalScore    = as.integer(MHRA_SIGNAL == "YES") +
                       as.integer(ROR_SIGNAL  == "YES") +
                       as.integer(IC_SIGNAL   == "YES") +
                       as.integer(EBGM_SIGNAL == "YES"),
      ACTIVE_SUBSTANCE = active_substance
    ) %>%
    arrange(desc(SignalScore), desc(PRR), desc(New))

  results
}

############################################################
# ESECUZIONE DEI PRODOTTI
############################################################

rmr_list <- list(

  # Fix: parentesi esplicite intorno a (eff_flag == 1 | str_detect(...,'BINOSTO'))
  # – senza parentesi, l'OR aveva precedenza su tutto il filtro sinistro
  RMR(
    paste0(
      "str_detect(GENERIC_NAME, 'ALENDRONATE SODIUM') & drug_name != '' & ",
      "!str_detect(drug_name, 'NOT COMPANY') & ",
      "(eff_flag == 1 | str_detect(drug_name, 'BINOSTO'))"
    ),
    "ALENDRONATE SODIUM (Form: Effervescent Tablet)"
  ),

  RMR(
    paste0(
      "str_detect(GENERIC_NAME, 'ALENDRONATE SODIUM') & drug_name != '' & ",
      "!str_detect(drug_name, 'NOT COMPANY') & ",
      "(is.na(eff_flag) | eff_flag != 1) & !str_detect(drug_name, 'BINOSTO')"
    ),
    "ALENDRONATE SODIUM (Form: All other forms)"
  ),

  RMR(drug_f("CLODRONATE DISODIUM"),                          "CLODRONATE DISODIUM"),
  RMR(drug_f("CLODRONATE DISODIUM, LIDOCAINE HYDROCHLORIDE"), "CLODRONATE DISODIUM/LIDOCAINE HYDROCHLORIDE"),
  RMR(drug_f("COLECALCIFEROL"),                               "COLECALCIFEROL"),
  RMR(drug_f("DIACEREIN"),                                    "DIACEREIN"),
  RMR(drug_f("ETORICOXIB"),                                   "ETORICOXIB"),
  RMR(drug_f("FLURBIPROFEN"),                                 "FLURBIPROFEN"),
  RMR(drug_f("GLIBENCLAMIDE"),                                "GLIBENCLAMIDE"),
  RMR(drug_f("LYOPHILIZED BACTERIAL LYSATES"),                "LYOPHILIZED BACTERIAL LYSATES"),
  RMR(drug_f("METFORMIN HYDROCHLORIDE"),                      "METFORMIN HYDROCHLORIDE"),
  RMR(drug_f("MOMETASONE FUROATE"),                           "MOMETASONE FUROATE"),
  RMR(drug_f("PARACETAMOL"),                                  "PARACETAMOL"),
  RMR(drug_f("PARACETAMOL, CODEINE PHOSPHATE"),               "PARACETAMOL/CODEINE PHOSPHATE"),
  RMR(drug_f("RRR-ALPHA-TOCOPHEROL"),                         "RRR-A-TOCOPHEROL"),
  RMR(drug_f("TACALCITOL MONOHYDRATE"),                       "TACALCITOL MONOHYDRATE"),
  RMR(drug_f("XYLOMETAZOLINE HYDROCHLORIDE"),                 "XYLOMETAZOLINE HYDROCHLORIDE")

)

############################################################
# RISULTATO FINALE
############################################################

final_results <- bind_rows(rmr_list)
View(final_results)

# write.csv(final_results, "RMR_Final_Results.csv", row.names = FALSE)
