/* ============================================================================
   sample_size_rct.sas
   -------------------
   Calcolo del campione per Studio Clinico Randomizzato (RCT) a due gruppi
   Outcome: variabile binaria (proporzione)

   Ipotesi statistiche:
     H0: p_trattamento = p_controllo  (nessuna differenza)
     H1: p_trattamento ≠ p_controllo  (differenza rilevabile)

   Parametri dello studio:
     - Proporzione attesa nel gruppo controllo (p0) : 30%
     - Potenza (1 - β)                             : 90%
     - Livello di significatività (α)              : 5% (two-sided)
     - Dropout / perdita al follow-up              : 10%
     - Test statistico                             : Chi-quadrato di Pearson

   Tested on SAS 9.4 / SAS Viya 4
   ============================================================================ */


/* ── 0. OPZIONI GLOBALI E MACRO VARIABILI ──────────────────────────────────── */

options nodate nonumber ls=120 ps=max;

%let P0      = 0.30;   /* proporzione outcome nel gruppo CONTROLLO            */
%let ALPHA   = 0.05;   /* livello di significatività (two-sided)              */
%let POWER   = 0.90;   /* potenza desiderata                                  */
%let DROPOUT = 0.10;   /* tasso di dropout / perdita al follow-up             */


/* ── 1. CALCOLO N BASE CON PROC POWER ─────────────────────────────────────── */
/*
   Analisi di sensibilità: vengono valutate differenze assolute rilevabili
   pari a +10%, +15%, +20%, +25%, +30% rispetto alla proporzione di controllo.
   (p_trattamento = 40%, 45%, 50%, 55%, 60%)
*/

proc power;
    title  "RCT Due Gruppi – Calcolo Sample Size per Proporzione Binaria";
    title2 "Controllo p0=30% | Potenza=90% | Alpha=5% | Two-sided | Test chi-quadrato";

    twosamplefreq test = pchi
        refproportion = &P0                        /* p controllo = 30%        */
        proportiondiff = 0.10 0.15 0.20 0.25 0.30 /* differenze da rilevare   */
        alpha          = &ALPHA
        power          = &POWER
        sides          = 2
        npergroup      = .;                        /* incognita: N per gruppo  */

    ods output output = _power_base;
run;
title;


/* ── 2. AGGIUSTAMENTO PER DROPOUT E CALCOLO N FINALE ─────────────────────── */

data WORK.SAMPLE_SIZE_FINAL;
    set _power_base;

    /* Proporzioni dei due gruppi */
    p_controllo   = &P0;
    p_trattamento = &P0 + ProportionDiff;

    /* Differenza assoluta (%) */
    diff_assoluta_pct = ProportionDiff * 100;

    /* N per gruppo prima della correzione dropout */
    n_per_gruppo_raw = NPerGroup;

    /* N per gruppo corretto per dropout (arrotondato per eccesso) */
    n_per_gruppo_adj = ceil(NPerGroup / (1 - &DROPOUT));

    /* N totale dello studio (due gruppi) */
    n_totale_adj = n_per_gruppo_adj * 2;

    label
        p_controllo       = "P Controllo"
        p_trattamento     = "P Trattamento"
        diff_assoluta_pct = "Differenza Assoluta (%)"
        NominalPower      = "Potenza Nominale"
        n_per_gruppo_raw  = "N/gruppo (senza dropout)"
        n_per_gruppo_adj  = "N/gruppo (dropout +10%)"
        n_totale_adj      = "N Totale (dropout +10%)";

    keep p_controllo p_trattamento diff_assoluta_pct NominalPower
         n_per_gruppo_raw n_per_gruppo_adj n_totale_adj;
run;


/* ── 3. REPORT TABELLARE FINALE ───────────────────────────────────────────── */

proc print data=WORK.SAMPLE_SIZE_FINAL noobs label;
    title  "Sample Size Finale – RCT Due Gruppi";
    title2 "Proporzione Controllo=30% | Potenza=90% | Alpha=5% | Dropout=10%";

    var p_controllo p_trattamento diff_assoluta_pct NominalPower
        n_per_gruppo_raw n_per_gruppo_adj n_totale_adj;

    format p_controllo p_trattamento   percent6.1
           diff_assoluta_pct           5.1
           NominalPower                percent6.1
           n_per_gruppo_raw
           n_per_gruppo_adj
           n_totale_adj                6.0;
run;
title;


/* ── 4. CURVA DI POTENZA ──────────────────────────────────────────────────── */
/*
   Grafico: N per gruppo (con dropout) vs Potenza, per varie differenze
   da rilevare. Utile per presentazioni e documenti di pianificazione.
*/

proc power plotonly;
    title  "Curva di Potenza – RCT Due Gruppi";
    title2 "Controllo p0=30% | Alpha=5% | Two-sided";

    twosamplefreq test = pchi
        refproportion = &P0
        proportiondiff = 0.10 0.15 0.20 0.25 0.30
        alpha          = &ALPHA
        power          = 0.70 to 0.95 by 0.01
        sides          = 2
        npergroup      = .;

    plot x=n step=1
         key=byfeature(pos=inset)
         yopts=(ref=&POWER refline=yes label="Power=90%")
         vary(linestyle by proportiondiff);
run;
title;


/* ── 5. VERIFICA: ONE-SAMPLE (opzionale) ──────────────────────────────────── */
/*
   Scenario alternativo: test di una sola proporzione vs valore di riferimento
   (es. testare se la proporzione osservata = 30% vs H0: p = p_null).
   Utile se lo studio ha un solo braccio attivo da confrontare con un benchmark.
*/

proc power;
    title  "One-Sample – Test Proporzione vs Valore di Riferimento (Benchmark)";
    title2 "H1: p=30% vs H0: p_null | Potenza=90% | Alpha=5%";

    onesamplefreq test = z
        proportion     = &P0             /* p osservato = 30%               */
        nullproportion = 0.20 0.25       /* valori null da testare          */
        alpha          = &ALPHA
        power          = &POWER
        sides          = 2
        ntotal         = .;

    ods output output = _power_1samp;
run;
title;

data WORK.ONE_SAMPLE_ADJ;
    set _power_1samp;

    n_adj     = ceil(NTotal / (1 - &DROPOUT));
    diff_null = (Proportion - NullProportion) * 100;

    label
        NullProportion = "P Nulla (H0)"
        Proportion     = "P Osservata (H1)"
        diff_null      = "Differenza da H0 (%)"
        NominalPower   = "Potenza Nominale"
        NTotal         = "N Totale (senza dropout)"
        n_adj          = "N Totale (dropout +10%)";

    keep NullProportion Proportion diff_null NominalPower NTotal n_adj;
run;

proc print data=WORK.ONE_SAMPLE_ADJ noobs label;
    title  "One-Sample: N con Correzione Dropout 10%";
    format NullProportion Proportion percent6.1
           diff_null                 5.1
           NominalPower              percent6.1
           NTotal n_adj              6.0;
run;
title;


/* ── 6. RIEPILOGO PARAMETRI IN LOG ───────────────────────────────────────── */

%put NOTE: =============================================================;
%put NOTE: PARAMETRI STUDIO;
%put NOTE:   Proporzione controllo (p0)  = &P0;
%put NOTE:   Livello alpha (two-sided)   = &ALPHA;
%put NOTE:   Potenza target              = &POWER;
%put NOTE:   Tasso di dropout            = &DROPOUT;
%put NOTE: =============================================================;


/* ── FINE PROGRAMMA ─────────────────────────────────────────────────────── */
