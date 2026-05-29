/* ============================================================================
   pharmacovigilance.sas
   ----------------------
   Pharmacovigilance analysis toolkit covering:
     1.  Sample data simulation (ICSR / adverse-event dataset)
     2.  Data cleaning & validation
     3.  Adverse Event (AE) frequency tables
     4.  Time-to-onset (TTO) analysis
     5.  Seriousness & outcome classification
     6.  Disproportionality analysis  (PRR, ROR, IC)
     7.  Signal detection summary
     8.  Automated case narrative generation

   Tested on SAS 9.4 / SAS Viya 4
   ============================================================================ */


/* ── 0. GLOBAL OPTIONS & MACRO VARIABLES ─────────────────────────────────── */

options nodate nonumber ls=120 ps=max;
%let PRODUCT   = DrugX;         /* product under surveillance            */
%let CUTDATE   = 31DEC2024;     /* data-lock date                        */
%let MINCOUNT  = 3;             /* min. case count for signal detection  */
%let PRR_THRESH= 2;             /* PRR threshold for flagging signals    */
%let CHI_THRESH= 4;             /* Chi-square threshold (Gamma Poisson)  */


/* ── 1. SIMULATE ICSR / AE INPUT DATA ────────────────────────────────────── */
/* In production, replace this step with a read from your safety database    */

data WORK.ICSR_RAW;
    length case_id        $12
           country        $3
           drug_name      $50
           ae_term        $80
           soc            $80   /* System Organ Class                   */
           seriousness    $3    /* YES / NO                             */
           outcome        $30
           reporter_type  $20
           sex            $1
           age_group      $20;

    infile datalines dsd missover;
    input case_id $ report_date :date9. onset_date :date9.
          country $ drug_name $ dose_mg ae_term $ soc $
          seriousness $ outcome $ reporter_type $ sex $ age $;
    format report_date onset_date date9.;
    age_group = age;   /* rename for clarity */
    drop age;

    /* time-to-onset in days */
    if onset_date ne . and report_date ne . then
        tto_days = report_date - onset_date;

datalines;
CASE0001,01JAN2024,28DEC2023,USA,DrugX,100,Nausea,Gastrointestinal disorders,NO,Recovered,Physician,F,Adult (18-65)
CASE0002,05JAN2024,03JAN2024,GBR,DrugX,200,Headache,Nervous system disorders,NO,Recovered,Consumer,M,Adult (18-65)
CASE0003,10JAN2024,07JAN2024,DEU,DrugX,100,Liver enzyme increased,Investigations,YES,Recovering,Physician,F,Elderly (>65)
CASE0004,12JAN2024,10JAN2024,FRA,DrugX,50,Rash,Skin and subcutaneous tissue disorders,NO,Recovered,Pharmacist,M,Adult (18-65)
CASE0005,15JAN2024,12JAN2024,USA,DrugX,200,Anaphylaxis,Immune system disorders,YES,Recovered,Physician,F,Adult (18-65)
CASE0006,18JAN2024,15JAN2024,CAN,DrugX,100,Dizziness,Nervous system disorders,NO,Recovered,Consumer,M,Elderly (>65)
CASE0007,20JAN2024,16JAN2024,USA,DrugX,200,Liver failure,Hepatobiliary disorders,YES,Fatal,Physician,F,Elderly (>65)
CASE0008,22JAN2024,20JAN2024,AUS,DrugX,100,Nausea,Gastrointestinal disorders,NO,Recovered,Consumer,M,Adult (18-65)
CASE0009,25JAN2024,22JAN2024,USA,DrugX,50,Fatigue,General disorders,NO,Recovered,Nurse,F,Adult (18-65)
CASE0010,28JAN2024,25JAN2024,ITA,DrugX,200,Thrombocytopenia,Blood and lymphatic system disorders,YES,Recovering,Physician,M,Elderly (>65)
CASE0011,01FEB2024,28JAN2024,USA,DrugX,100,Nausea,Gastrointestinal disorders,NO,Recovered,Consumer,F,Adult (18-65)
CASE0012,03FEB2024,31JAN2024,GBR,DrugX,200,Anaphylaxis,Immune system disorders,YES,Recovered,Physician,M,Adult (18-65)
CASE0013,07FEB2024,04FEB2024,USA,DrugX,100,Headache,Nervous system disorders,NO,Recovered,Consumer,F,Pediatric (<18)
CASE0014,10FEB2024,07FEB2024,DEU,DrugX,50,Liver enzyme increased,Investigations,YES,Recovering,Physician,M,Adult (18-65)
CASE0015,14FEB2024,10FEB2024,USA,DrugX,200,Rash,Skin and subcutaneous tissue disorders,NO,Recovered,Consumer,F,Adult (18-65)
;
run;


/* ── 2. DATA CLEANING & VALIDATION ──────────────────────────────────────── */

/* 2a. Flag missing critical fields */
data WORK.ICSR_CLEAN WORK.ICSR_ERRORS;
    set WORK.ICSR_RAW;

    length error_flags $200;
    error_flags = '';

    if case_id     = '' then error_flags = catx('|', error_flags, 'MISSING_CASE_ID');
    if ae_term     = '' then error_flags = catx('|', error_flags, 'MISSING_AE');
    if drug_name   = '' then error_flags = catx('|', error_flags, 'MISSING_DRUG');
    if seriousness = '' then error_flags = catx('|', error_flags, 'MISSING_SERIOUSNESS');
    if report_date = .  then error_flags = catx('|', error_flags, 'MISSING_REPORT_DATE');
    if report_date > input("&CUTDATE", date9.)
                        then error_flags = catx('|', error_flags, 'POST_CUTOFF');
    if tto_days < 0     then error_flags = catx('|', error_flags, 'NEGATIVE_TTO');

    if error_flags ne '' then output WORK.ICSR_ERRORS;
    else                      output WORK.ICSR_CLEAN;
run;

/* 2b. Validation report */
proc sql;
    title "Data Validation Summary";
    select count(*) as total_raw,
           (select count(*) from WORK.ICSR_CLEAN)  as clean,
           (select count(*) from WORK.ICSR_ERRORS) as errors
    from WORK.ICSR_RAW;
quit;
title;

proc print data=WORK.ICSR_ERRORS noobs label;
    title "Records Failing Validation";
    var case_id ae_term error_flags;
run;
title;


/* ── 3. ADVERSE EVENT FREQUENCY TABLES ──────────────────────────────────── */

/* 3a. AE frequency by preferred term */
proc freq data=WORK.ICSR_CLEAN order=freq;
    title "AE Frequency by Preferred Term";
    tables ae_term / nocum nopercent;
run;
title;

/* 3b. AE frequency by System Organ Class */
proc freq data=WORK.ICSR_CLEAN order=freq;
    title "AE Frequency by System Organ Class (SOC)";
    tables soc / nocum nopercent;
run;
title;

/* 3c. Seriousness breakdown */
proc freq data=WORK.ICSR_CLEAN;
    title "Seriousness Distribution";
    tables seriousness * outcome / norow nopct;
run;
title;

/* 3d. Demographics cross-tab */
proc freq data=WORK.ICSR_CLEAN;
    title "Cases by Sex and Age Group";
    tables sex * age_group / norow nopct;
run;
title;

/* 3e. Country distribution */
proc freq data=WORK.ICSR_CLEAN order=freq;
    title "Case Distribution by Country";
    tables country / nocum nopercent;
run;
title;


/* ── 4. TIME-TO-ONSET ANALYSIS ──────────────────────────────────────────── */

proc means data=WORK.ICSR_CLEAN n mean median min max std
           maxdec=1;
    title "Time-to-Onset (days) – Overall";
    where tto_days ne .;
    var tto_days;
run;
title;

proc means data=WORK.ICSR_CLEAN n mean median min max
           maxdec=1;
    title "Time-to-Onset (days) by SOC";
    where tto_days ne .;
    class soc;
    var tto_days;
run;
title;

/* Histogram of TTO */
proc sgplot data=WORK.ICSR_CLEAN;
    title "Time-to-Onset Distribution";
    where tto_days ne .;
    histogram tto_days / binwidth=2 fillattrs=(color=steelblue);
    density   tto_days / type=kernel lineattrs=(color=red);
    xaxis label="Days from Drug Start to AE Onset";
    yaxis label="Number of Cases";
run;
title;


/* ── 5. SERIOUSNESS & OUTCOME CLASSIFICATION ────────────────────────────── */

/* Summarise serious vs non-serious by AE term */
proc tabulate data=WORK.ICSR_CLEAN;
    title "Serious vs Non-Serious Cases by AE Term";
    class ae_term seriousness;
    table ae_term, seriousness*(N) / rts=50;
run;
title;

/* Fatal cases */
proc print data=WORK.ICSR_CLEAN noobs label;
    title "Fatal Cases";
    where upcase(outcome) = 'FATAL';
    var case_id report_date country drug_name dose_mg ae_term soc age_group sex;
run;
title;


/* ── 6. DISPROPORTIONALITY ANALYSIS ─────────────────────────────────────── */
/*
   Metrics computed per AE term:
     a11 = cases with PRODUCT + AE
     a12 = cases with PRODUCT + other AE
     a21 = cases with other product + AE      (requires background data)
     a22 = cases with other product + other AE

   For demonstration we construct a simplified 2×2 table using simulated
   background counts.  In production, populate N_BACKGROUND from your
   full pharmacovigilance database (e.g., WHO-VigiBase, FAERS).
*/

%let N_TOTAL_DB   = 500000;   /* total reports in background database  */
%let N_PRODUCT_DB = 15;       /* total reports for DrugX in database   */

proc sql;
    create table WORK.DISP as
    select
        ae_term,
        soc,
        count(*)                as a11,          /* product + AE          */
        &N_PRODUCT_DB - count(*) as a12,         /* product + other AE    */
        /* simulated background (in real use: query your DB) */
        max(2, int(count(*) * 15))  as a21,      /* other drug + same AE  */
        &N_TOTAL_DB - &N_PRODUCT_DB
          - max(2, int(count(*)*15)) as a22      /* other drug + other AE */
    from WORK.ICSR_CLEAN
    group by ae_term, soc
    having count(*) >= &MINCOUNT
    order by a11 desc;
quit;

/* Compute PRR, ROR, and IC (Information Component) */
data WORK.SIGNALS;
    set WORK.DISP;

    /* Proportional Reporting Ratio */
    if (a11 + a12) > 0 and (a21 + a22) > 0 then do;
        prr      = (a11 / (a11 + a12)) / (a21 / (a21 + a22));
        log_prr  = log(prr);
        se_log_prr = sqrt(1/a11 - 1/(a11+a12) + 1/a21 - 1/(a21+a22));
        prr_lower95 = exp(log_prr - 1.96 * se_log_prr);
        prr_upper95 = exp(log_prr + 1.96 * se_log_prr);
    end;

    /* Reporting Odds Ratio */
    if a12 > 0 and a21 > 0 then do;
        ror      = (a11 * a22) / (a12 * a21);
        se_log_ror = sqrt(1/a11 + 1/a12 + 1/a21 + 1/a22);
        ror_lower95 = exp(log(ror) - 1.96 * se_log_ror);
        ror_upper95 = exp(log(ror) + 1.96 * se_log_ror);
    end;

    /* Information Component (Bayesian; simplified observed/expected) */
    n_total = a11 + a12 + a21 + a22;
    expected = ((a11 + a12) * (a11 + a21)) / n_total;
    if expected > 0 then
        ic = log2(a11 / expected);

    /* Chi-square statistic */
    if n_total > 0 then
        chi2 = n_total * ((a11*a22 - a12*a21)**2) /
               ((a11+a12)*(a21+a22)*(a11+a21)*(a12+a22));

    /* Signal flag */
    signal_flag = (prr >= &PRR_THRESH and a11 >= &MINCOUNT and
                   prr_lower95 >= 1 and chi2 >= &CHI_THRESH);

    label
        ae_term     = "Adverse Event Term"
        a11         = "Product+AE (n)"
        prr         = "PRR"
        prr_lower95 = "PRR 95% CI Lower"
        prr_upper95 = "PRR 95% CI Upper"
        ror         = "ROR"
        ror_lower95 = "ROR 95% CI Lower"
        ror_upper95 = "ROR 95% CI Upper"
        ic          = "IC (log2)"
        chi2        = "Chi-Square"
        signal_flag = "Signal Detected";
run;


/* ── 7. SIGNAL DETECTION SUMMARY ────────────────────────────────────────── */

proc print data=WORK.SIGNALS noobs label;
    title "Disproportionality Analysis – All AE Terms (n >= &MINCOUNT)";
    var ae_term soc a11 prr prr_lower95 prr_upper95
        ror ror_lower95 ror_upper95 ic chi2 signal_flag;
    format prr prr_lower95 prr_upper95
           ror ror_lower95 ror_upper95
           ic chi2  7.2;
run;
title;

proc print data=WORK.SIGNALS(where=(signal_flag=1)) noobs label;
    title "*** POTENTIAL SIGNALS (PRR >= &PRR_THRESH, Chi2 >= &CHI_THRESH, n >= &MINCOUNT) ***";
    var ae_term soc a11 prr prr_lower95 prr_upper95 ic chi2;
    format prr prr_lower95 prr_upper95 ic chi2  7.2;
run;
title;

/* Forest plot of PRR */
proc sgplot data=WORK.SIGNALS;
    title "PRR Forest Plot by AE Term";
    scatter y=ae_term x=prr /
        xerrorlower=prr_lower95 xerrorupper=prr_upper95
        markerattrs=(symbol=circlefilled size=8);
    refline 1 / axis=x lineattrs=(pattern=dash color=red);
    refline &PRR_THRESH / axis=x lineattrs=(pattern=shortdash color=orange)
            label="PRR threshold (&PRR_THRESH)";
    xaxis label="Proportional Reporting Ratio (95% CI)" min=0;
    yaxis label=" " discreteorder=data;
run;
title;


/* ── 8. AUTOMATED CASE NARRATIVE GENERATION ─────────────────────────────── */

%macro generate_narrative(dsn=WORK.ICSR_CLEAN, outfile=narratives.txt);
    filename narr "&outfile";
    data _null_;
        set &dsn;
        file narr;

        /* Format dates for readability */
        rd = put(report_date, worddate18.);
        od = put(onset_date,  worddate18.);

        /* Build narrative string */
        narrative = catx(' ',
            'Case', case_id, '–',
            'A', age_group, strip(sex='F' ? 'female' : 'male'),
            'patient from', country,
            'treated with', drug_name, '(' || strip(put(dose_mg,8.)) || ' mg)',
            'reported', ae_term, 'on', strip(od) || '.',
            'The case was classified as',
            (seriousness='YES' ? 'SERIOUS' : 'non-serious') || '.',
            'Outcome:', strip(outcome) || '.',
            'Reporter type:', strip(reporter_type) || '.',
            'Report received:', strip(rd) || '.'
        );
        if tto_days ne . then
            narrative = catx(' ', narrative,
                'Time to onset:', strip(put(tto_days, 8.)), 'days.');

        put narrative;
        put '---';
    run;
    filename narr clear;
    %put NOTE: Narratives written to &outfile;
%mend generate_narrative;

%generate_narrative(dsn=WORK.ICSR_CLEAN, outfile=narratives.txt);


/* ── 9. EXPORT SIGNAL TABLE TO EXCEL ───────────────────────────────────── */

ods excel file="pharmacovigilance_signals.xlsx"
    options(sheet_name="Signals" embedded_titles="yes");

proc print data=WORK.SIGNALS noobs label;
    title "Pharmacovigilance Signal Detection – &PRODUCT – Data lock: &CUTDATE";
    var ae_term soc a11 prr prr_lower95 prr_upper95
        ror ic chi2 signal_flag;
    format prr prr_lower95 prr_upper95 ror ic chi2  7.2;
run;

ods excel close;


/* ── 10. NERIDRONATE OFF-LABEL USE EXTRACTION (ORACLE PASS-THROUGH) ─────── */
/*  Requires macro vars: &ora_user, &ora_pass, &ora_path                       */
/*  CLOB fix: cp.generic_name accessed via dbms_lob.substr(…,200)              */

proc sql;
    connect to oracle (user=&ora_user password=&ora_pass path=&ora_path);

    create table WORK.NERIDRONATE_OFFLABEL as
    select * from connection to oracle
    (
        select distinct
            cm.case_id,
            cm.case_num,
            ce.seq_num,

            lg.gender,
            cpi.pat_age,
            cpi.ind_pref_term,
            lag.group_name                        as age_group,

            /* CLOB fix */
            dbms_lob.substr(cp.generic_name, 200) as generic_name,

            /* Init receipt date filter flag */
            case
                when cm.init_rept_date <= DATE '2025-07-04'
                then 'ok'
                else 'not ok'
            end                                   as init_receipt_date,

            /* Off-label event term (null for non-off-label rows) */
            case
                when lower(ce.pref_term) = 'off label use'
                then ce.pref_term
                else null
            end                                   as PT,

            ce.body_sys                           as soc

        from case_master cm

            inner join case_pat_info cpi
                on  cm.case_id      = cpi.case_id
                and cpi.deleted     is null

            left join lm_age_groups lag
                on  lag.age_group_id = cpi.age_group_id

            left join lm_gender lg
                on  lg.gender_id    = cpi.gender_id

            left join case_assess ca
                on  cm.case_id      = ca.case_id
                and ca.deleted      is null

            inner join case_product cp
                on  cm.case_id      = cp.case_id
                and cp.deleted      is null

            inner join case_event ce
                on  cm.case_id      = ce.case_id
                and ce.deleted      is null

        where cm.deleted is null
          and cm.state_id <> 1

          /* Off-label use filter */
          and lower(ce.pref_term) = 'off label use'

          /* Drug filter – CLOB fix */
          and upper(dbms_lob.substr(cp.generic_name, 200)) like '%NERIDRONATE SODIUM%'
          and upper(dbms_lob.substr(cp.generic_name, 200)) not like '%NOT COMPANY%'

        order by cm.case_id, ce.seq_num
    );

    disconnect from oracle;
quit;

/* Persist to network library */
libname mylib "\\Vivaldi\ShBiostatistics\FV\Neridr_pediatrico\26 Maggio 2026";

data mylib.NERIDRONATE_OFFLABEL;
    set WORK.NERIDRONATE_OFFLABEL;
run;

proc export data=mylib.NERIDRONATE_OFFLABEL
    outfile="\\Vivaldi\ShBiostatistics\FV\Neridr_pediatrico\26 Maggio 2026\NERIDRONATE_OFFLABEL.xlsx"
    dbms=xlsx
    replace;
run;

libname mylib clear;


/* ── END OF PROGRAM ─────────────────────────────────────────────────────── */
