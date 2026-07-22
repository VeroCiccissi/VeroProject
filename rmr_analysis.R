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
############################################################

interval     <- cumulative %>% filter(recdate >= datefrom)
lastinterval <- cumulative %>% filter(recdate >= lastdatefrom, recdate <= lastdateto)

############################################################
# FUNZIONE RMR
# FIX ROOT CAUSE: sostituito parse_expr(drug_filter) con enquo(drug_filter).
# parse_expr() richiede una stringa; passare GENERIC_NAME senza virgolette
# causava la valutazione immediata nell'ambiente globale → errore
# "oggetto 'GENERIC_NAME' non trovato".
# enquo() cattura l'espressione pigramente e la risolve dentro filter()
# dove GENERIC_NAME esiste come colonna del dataframe.
############################################################

RMR <- function(drug_filter, active_substance) {

  drug_quo          <- enquo(drug_filter)          # lazy capture – nessuna valutazione immediata
  cumulative_drug   <- cumulative   %>% filter(!!drug_quo)
  interval_drug     <- interval     %>% filter(!!drug_quo)
  lastinterval_drug <- lastinterval %>% filter(!!drug_quo)

  n_cum  <- n_distinct(cumulative_drug$CASE_ID)
  n_intv <- n_distinct(interval_drug$CASE_ID)

  if (n_cum == 0 || n_intv == 0) {
    return(tibble(ACTIVE_SUBSTANCE = active_substance, RESULT = "NO DATA FOUND"))
  }

  # ── Conteggi tabella 2x2 ─────────────────────────────────────────────────
  N   <- n_distinct(cumulative$CASE_ID)
  ApB <- n_cum

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

  # ── Metriche (vettorizzate – niente rowwise) ──────────────────────────────
  metrics <- a %>%
    left_join(apc, by = c("body_system", "preferred_event_term")) %>%
    mutate(
      B = ApB - A,
      C = ApC - A,
      D = (N - ApB) - C,

      Observed = A,
      Expected = ((A + B) * (A + C)) / (A + B + C + D),
      OE_Ratio = Observed / Expected,

      PRR     = if_else(C > 0, (A / (A + B)) / (C / (C + D)),        NA_real_),
      PRR_SE  = if_else(C > 0, sqrt(1/A - 1/(A+B) + 1/C - 1/(C+D)), NA_real_),
      PRR_LCL = exp(log(PRR) - 1.96 * PRR_SE),
      PRR_UCL = exp(log(PRR) + 1.96 * PRR_SE),

      ROR     = if_else(B > 0 & C > 0 & D > 0, (A * D) / (B * C),            NA_real_),
      ROR_SE  = if_else(B > 0 & C > 0 & D > 0, sqrt(1/A + 1/B + 1/C + 1/D), NA_real_),
      ROR_LCL = exp(log(ROR) - 1.96 * ROR_SE),
      ROR_UCL = exp(log(ROR) + 1.96 * ROR_SE),

      CHISQ = ((A+B+C+D) * (A*D - B*C)^2) / ((A+B) * (C+D) * (A+C) * (B+D)),

      IC    = log2((Observed + 0.5) / (Expected + 0.5)),
      IC_SE = sqrt(1 / (Observed + 0.5)),
      IC025 = IC - 1.96 * IC_SE,
      IC975 = IC + 1.96 * IC_SE,

      EBGM  = (Observed + 0.5) / (Expected + 0.5)
    ) %>%
    mutate(
      FISHER_P = pmap_dbl(
        list(A, B, C, D),
        function(a, b, c, d) {
          tryCatch(
            fisher.test(matrix(c(a, b, c, d), nrow = 2))$p.value,
            error = function(e) NA_real_
          )
        }
      ),
      MHRA_SIGNAL = if_else(!is.na(PRR) & A >= 3 & PRR >= 2 & CHISQ >= 4, "YES", "NO"),
      ROR_SIGNAL  = if_else(!is.na(ROR_LCL) & ROR_LCL > 1,                "YES", "NO"),
      IC_SIGNAL   = if_else(!is.na(IC025)   & IC025   > 0,                "YES", "NO"),
      EBGM_SIGNAL = if_else(!is.na(EBGM)    & EBGM   >= 2,                "YES", "NO")
    ) %>%
    rename(Total = A)

  # ── Risultati finali ──────────────────────────────────────────────────────
  results <- metrics %>%
    left_join(intvcount_PT,  by = c("body_system", "preferred_event_term")) %>%
    left_join(lastcount_PT,  by = c("body_system", "preferred_event_term")) %>%
    left_join(fatal_PT,      by = c("body_system", "preferred_event_term")) %>%
    left_join(intvfatal_PT,  by = c("body_system", "preferred_event_term")) %>%
    replace_na(list(New = 0L, Last = 0L, TotalFatal = 0L, NewFatal = 0L)) %>%
    filter(New > 0) %>%
    mutate(
      GrowthRatio      = round(New / (Last + 1), 2),
      QoQ_Percent      = round(100 * (New - Last) / pmax(Last, 1), 2),
      FatalityRatio    = round(100 * TotalFatal / Total, 2),
      SignalScore      = as.integer(MHRA_SIGNAL == "YES") +
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
# Le espressioni sono passate SENZA virgolette – enquo() le cattura
# correttamente e le risolve nel contesto del dataframe.
############################################################

rmr_list <- list(

  RMR(
    str_detect(GENERIC_NAME, 'ALENDRONATE SODIUM') & drug_name != '' &
      !str_detect(drug_name, 'NOT COMPANY') &
      (eff_flag == 1 | str_detect(drug_name, 'BINOSTO')),
    "ALENDRONATE SODIUM (Form: Effervescent Tablet)"
  ),

  RMR(
    str_detect(GENERIC_NAME, 'ALENDRONATE SODIUM') & drug_name != '' &
      !str_detect(drug_name, 'NOT COMPANY') &
      (is.na(eff_flag) | eff_flag != 1) & !str_detect(drug_name, 'BINOSTO'),
    "ALENDRONATE SODIUM (Form: All other forms)"
  ),

  RMR(GENERIC_NAME == 'CLODRONATE DISODIUM'                        & drug_name != '' & !str_detect(drug_name, 'NOT COMPANY'), "CLODRONATE DISODIUM"),
  RMR(GENERIC_NAME == 'CLODRONATE DISODIUM, LIDOCAINE HYDROCHLORIDE' & drug_name != '' & !str_detect(drug_name, 'NOT COMPANY'), "CLODRONATE DISODIUM/LIDOCAINE HYDROCHLORIDE"),
  RMR(GENERIC_NAME == 'COLECALCIFEROL'                              & drug_name != '' & !str_detect(drug_name, 'NOT COMPANY'), "COLECALCIFEROL"),
  RMR(GENERIC_NAME == 'DIACEREIN'                                   & drug_name != '' & !str_detect(drug_name, 'NOT COMPANY'), "DIACEREIN"),
  RMR(GENERIC_NAME == 'ETORICOXIB'                                  & drug_name != '' & !str_detect(drug_name, 'NOT COMPANY'), "ETORICOXIB"),
  RMR(GENERIC_NAME == 'FLURBIPROFEN'                                & drug_name != '' & !str_detect(drug_name, 'NOT COMPANY'), "FLURBIPROFEN"),
  RMR(GENERIC_NAME == 'GLIBENCLAMIDE'                               & drug_name != '' & !str_detect(drug_name, 'NOT COMPANY'), "GLIBENCLAMIDE"),
  RMR(GENERIC_NAME == 'LYOPHILIZED BACTERIAL LYSATES'               & drug_name != '' & !str_detect(drug_name, 'NOT COMPANY'), "LYOPHILIZED BACTERIAL LYSATES"),
  RMR(GENERIC_NAME == 'METFORMIN HYDROCHLORIDE'                     & drug_name != '' & !str_detect(drug_name, 'NOT COMPANY'), "METFORMIN HYDROCHLORIDE"),
  RMR(GENERIC_NAME == 'MOMETASONE FUROATE'                         & drug_name != '' & !str_detect(drug_name, 'NOT COMPANY'), "MOMETASONE FUROATE"),
  RMR(GENERIC_NAME == 'PARACETAMOL'                                 & drug_name != '' & !str_detect(drug_name, 'NOT COMPANY'), "PARACETAMOL"),
  RMR(GENERIC_NAME == 'PARACETAMOL, CODEINE PHOSPHATE'              & drug_name != '' & !str_detect(drug_name, 'NOT COMPANY'), "PARACETAMOL/CODEINE PHOSPHATE"),
  RMR(GENERIC_NAME == 'RRR-ALPHA-TOCOPHEROL'                        & drug_name != '' & !str_detect(drug_name, 'NOT COMPANY'), "RRR-A-TOCOPHEROL"),
  RMR(GENERIC_NAME == 'TACALCITOL MONOHYDRATE'                      & drug_name != '' & !str_detect(drug_name, 'NOT COMPANY'), "TACALCITOL MONOHYDRATE"),
  RMR(GENERIC_NAME == 'XYLOMETAZOLINE HYDROCHLORIDE'                & drug_name != '' & !str_detect(drug_name, 'NOT COMPANY'), "XYLOMETAZOLINE HYDROCHLORIDE")

)

############################################################
# RISULTATO FINALE
############################################################

final_results <- bind_rows(rmr_list)
View(final_results)

# write.csv(final_results, "RMR_Final_Results.csv", row.names = FALSE)
