/* ============================================================================
   sample_size_rct.sas
   -------------------
   Calcolo del campione per Studio Clinico Randomizzato (RCT) a due gruppi
   Outcome: variabile binaria (proporzione)

   DEFINIZIONE ENDPOINT PRIMARIO
   ─────────────────────────────
   Proporzione di pazienti con ≥ 5 side effects (SE)
   nei 90 giorni successivi all'ultima somministrazione del farmaco.

   BASI STORICHE (studio precedente)
   ──────────────────────────────────
   I valori di riferimento derivano da uno studio precedente che ha osservato:
     - Gruppo A (controllo) : 44%  dei pazienti con ≥5 SE in 90 gg
     - Gruppo B (trattamento): 21%  dei pazienti con ≥5 SE in 90 gg
     - Proporzione ponderata (pooled reference): 30%

   IPOTESI STATISTICHE
   ────────────────────
     H0 : p_trattamento = p_controllo   (nessuna differenza tra i gruppi)
     H1 : p_trattamento ≠ p_controllo   (differenza rilevabile, two-sided)

   PARAMETRI DELLO STUDIO
   ───────────────────────
     - Proporzione gruppo controllo  (p_ctrl)   : 44%   [dati storici]
     - Proporzione gruppo trattamento (p_treat)  : 21%   [dati storici]
     - Proporzione ponderata pooled  (p_pool)   : 30%   [riferimento]
     - Potenza (1 − β)                          : 90%
     - Livello di significatività (α)           : 5%  (two-sided)
     - Dropout / perdita al follow-up           : 10%
     - Test statistico                          : Chi-quadrato di Pearson (PCHI)
     - Finestra di osservazione                 : 90 giorni post ultima dose
     - Soglia evento                            : ≥ 5 side effects

   Tested on SAS 9.4 / SAS Viya 4
   ============================================================================ */


/* ── 0. OPZIONI GLOBALI E MACRO VARIABILI ──────────────────────────────────── */

options nodate nonumber ls=120 ps=max;

/* Proporzioni osservate nello studio storico */
%let P_CTRL    = 0.44;   /* gruppo controllo  – 44% dei pz con ≥5 SE in 90gg  */
%let P_TREAT   = 0.21;   /* gruppo trattamento – 21% dei pz con ≥5 SE in 90gg  */
%let P_POOL    = 0.30;   /* proporzione ponderata pooled (riferimento)          */

/* Parametri dello studio */
%let ALPHA     = 0.05;   /* livello di significatività (two-sided)              */
%let POWER     = 0.90;   /* potenza desiderata                                  */
%let DROPOUT   = 0.10;   /* tasso di dropout / perdita al follow-up             */

/* Definizione endpoint */
%let MIN_SE    = 5;      /* numero minimo di SE per definire l'evento           */
%let FOLLOWUP  = 90;     /* giorni di osservazione post ultima somministrazione */


/* ── 1. VERIFICA PROPORZIONE PONDERATA (dati storici) ─────────────────────── */
/*
   Calcola il peso implicito dei due bracci storici che produce p_pool = 30%.
   Formula: p_pool = w * p_ctrl + (1 - w) * p_treat
            w = (p_pool - p_treat) / (p_ctrl - p_treat)
*/

data _null_;
    p_ctrl  = &P_CTRL;
    p_treat = &P_TREAT;
    p_pool  = &P_POOL;

    /* Peso implicito del gruppo controllo nello studio storico */
    w_ctrl  = (p_pool - p_treat) / (p_ctrl - p_treat);
    w_treat = 1 - w_ctrl;

    /* Verifica aritmetica */
    p_pool_check = w_ctrl * p_ctrl + w_treat * p_treat;

    put '─────────────────────────────────────────────────────';
    put 'VERIFICA PROPORZIONE PONDERATA (studio storico)';
    put '─────────────────────────────────────────────────────';
    put 'Gruppo controllo  : p=' p_ctrl  percent6.1 '  peso=' w_ctrl  percent6.1;
    put 'Gruppo trattamento: p=' p_treat percent6.1 '  peso=' w_treat percent6.1;
    put 'Proporzione pooled calcolata : ' p_pool_check percent6.1;
    put 'Proporzione pooled dichiarata: ' p_pool       percent6.1;
    put '─────────────────────────────────────────────────────';
run;


/* ── 2. ANALISI PRIMARIA – PROPORZIONI STORICHE (44% vs 21%) ──────────────── */
/*
   Calcolo del campione con le proporzioni osservate nello studio precedente.
   Rappresenta la stima principale del fabbisogno di campione.
*/

proc power;
    title  "RCT Due Gruppi – Analisi Primaria (dati storici)";
    title2 "Endpoint: ≥&MIN_SE SE nei &FOLLOWUP gg post ultima dose";
    title3 "Controllo=44% vs Trattamento=21% | Power=90% | Alpha=5% | Two-sided";

    twosamplefreq test = pchi
        groupproportions = (&P_CTRL &P_TREAT)   /* 0.44 vs 0.21               */
        alpha            = &ALPHA
        power            = &POWER
        sides            = 2
        npergroup        = .;                   /* incognita: N per gruppo     */

    ods output output = _power_primary;
run;
title;


/* ── 3. AGGIUSTAMENTO DROPOUT – ANALISI PRIMARIA ─────────────────────────── */

data WORK.SS_PRIMARY;
    set _power_primary;

    p_controllo      = Proportion1;             /* 44%                         */
    p_trattamento    = Proportion2;             /* 21%                         */
    p_pooled_ref     = &P_POOL;                 /* 30% – riferimento ponderato */
    diff_assoluta    = (Proportion1 - Proportion2) * 100; /* differenza in pp  */

    n_per_gruppo_raw = NPerGroup;
    n_per_gruppo_adj = ceil(NPerGroup / (1 - &DROPOUT));  /* correzione 10%    */
    n_totale_adj     = n_per_gruppo_adj * 2;

    label
        p_controllo      = "P Controllo (storico)"
        p_trattamento    = "P Trattamento (storico)"
        p_pooled_ref     = "P Ponderata Pooled"
        diff_assoluta    = "Differenza Assoluta (pp)"
        NominalPower     = "Potenza Nominale"
        n_per_gruppo_raw = "N/gruppo (senza dropout)"
        n_per_gruppo_adj = "N/gruppo (dropout +10%)"
        n_totale_adj     = "N Totale (dropout +10%)";

    keep p_controllo p_trattamento p_pooled_ref diff_assoluta NominalPower
         n_per_gruppo_raw n_per_gruppo_adj n_totale_adj;
run;

proc print data=WORK.SS_PRIMARY noobs label;
    title  "ANALISI PRIMARIA – Sample Size Finale";
    title2 "Endpoint: ≥&MIN_SE SE in &FOLLOWUP gg | Controllo=44% vs Trattamento=21%";
    title3 "Proporzione Ponderata Pooled (riferimento): &P_POOL | Dropout=10%";

    var p_controllo p_trattamento p_pooled_ref diff_assoluta NominalPower
        n_per_gruppo_raw n_per_gruppo_adj n_totale_adj;

    format p_controllo p_trattamento p_pooled_ref   percent6.1
           diff_assoluta                             6.1
           NominalPower                              percent6.1
           n_per_gruppo_raw n_per_gruppo_adj
           n_totale_adj                              6.0;
run;
title;


/* ── 4. ANALISI DI SENSIBILITÀ – VARIAZIONE PROPORZIONE TRATTAMENTO ──────── */
/*
   Si mantiene il gruppo controllo a p=44% (dato storico).
   Si varia la proporzione del gruppo trattamento da 0.15 a 0.35
   per valutare la robustezza del calcolo del campione a diverse
   ipotesi sull'efficacia del trattamento.

   Proporzione di riferimento trattamento (studio storico): 21%
   Differenze corrispondenti rispetto al controllo: da -9 pp a -29 pp
*/

proc power;
    title  "RCT Due Gruppi – Analisi di Sensibilità";
    title2 "Controllo fisso=44% | Variazione Proporzione Trattamento | Power=90% | Alpha=5%";

    twosamplefreq test = pchi
        refproportion  = &P_CTRL                       /* controllo = 44%      */
        proportiondiff = -0.09 -0.14 -0.19 -0.23 -0.29 /* da -9pp a -29pp     */
        /* corrispondenti a p_tratt = 0.35, 0.30, 0.25, 0.21*, 0.15           */
        /* * = valore storico principale                                        */
        alpha          = &ALPHA
        power          = &POWER
        sides          = 2
        npergroup      = .;

    ods output output = _power_sens;
run;
title;


/* ── 5. AGGIUSTAMENTO DROPOUT – ANALISI DI SENSIBILITÀ ───────────────────── */

data WORK.SS_SENSITIVITY;
    set _power_sens;

    p_controllo      = &P_CTRL;
    p_trattamento    = &P_CTRL + ProportionDiff;         /* p_ctrl + diff neg  */
    diff_assoluta    = abs(ProportionDiff) * 100;
    storico_flag     = (abs(ProportionDiff - (-0.23)) < 0.001); /* marca 21%  */

    n_per_gruppo_raw = NPerGroup;
    n_per_gruppo_adj = ceil(NPerGroup / (1 - &DROPOUT));
    n_totale_adj     = n_per_gruppo_adj * 2;

    label
        p_controllo      = "P Controllo"
        p_trattamento    = "P Trattamento"
        diff_assoluta    = "Riduzione Assoluta (pp)"
        storico_flag     = "Valore Storico"
        NominalPower     = "Potenza Nominale"
        n_per_gruppo_raw = "N/gruppo (senza dropout)"
        n_per_gruppo_adj = "N/gruppo (dropout +10%)"
        n_totale_adj     = "N Totale (dropout +10%)";

    keep p_controllo p_trattamento diff_assoluta storico_flag NominalPower
         n_per_gruppo_raw n_per_gruppo_adj n_totale_adj;
run;

proc print data=WORK.SS_SENSITIVITY noobs label;
    title  "ANALISI DI SENSIBILITÀ – Sample Size per Vari Scenari";
    title2 "Controllo=44% fisso | Endpoint: ≥&MIN_SE SE in &FOLLOWUP gg | Dropout=10%";
    title3 "(*) = scenario storico principale: Trattamento=21%";

    var p_controllo p_trattamento diff_assoluta storico_flag NominalPower
        n_per_gruppo_raw n_per_gruppo_adj n_totale_adj;

    format p_controllo p_trattamento   percent6.1
           diff_assoluta               6.1
           storico_flag                1.0
           NominalPower                percent6.1
           n_per_gruppo_raw
           n_per_gruppo_adj
           n_totale_adj                6.0;
run;
title;


/* ── 6. CURVA DI POTENZA ──────────────────────────────────────────────────── */
/*
   Grafico: N per gruppo vs Potenza per i diversi scenari di sensibilità.
   La linea verticale a Power=90% mostra il target dello studio.
*/

proc power plotonly;
    title  "Curva di Potenza – RCT Due Gruppi";
    title2 "Controllo=44% | Endpoint: ≥&MIN_SE SE in &FOLLOWUP gg | Alpha=5%";

    twosamplefreq test = pchi
        refproportion  = &P_CTRL
        proportiondiff = -0.09 -0.14 -0.19 -0.23 -0.29
        alpha          = &ALPHA
        power          = 0.70 to 0.95 by 0.01
        sides          = 2
        npergroup      = .;

    plot x     = n
         step  = 1
         key   = byfeature(pos=inset)
         yopts = (ref=&POWER refline=yes label="Power 90%")
         vary(linestyle by proportiondiff);
run;
title;


/* ── 7. TABELLA RIASSUNTIVA PARAMETRI FINALI ──────────────────────────────── */

data WORK.DESIGN_SUMMARY;
    length parametro $50 valore $30;
    infile datalines dsd;
    input parametro $ valore $;
datalines;
Endpoint primario,Proporzione pz con >=5 SE in 90gg post ultima dose
Proporzione gruppo controllo (storico),44%
Proporzione gruppo trattamento (storico),21%
Proporzione ponderata pooled (riferimento),30%
Differenza assoluta da rilevare,23 punti percentuali
Livello di significativita alpha,5% (two-sided)
Potenza (1-beta),90%
Test statistico,Chi-quadrato di Pearson (PCHI)
Tasso dropout,10%
Finestra di osservazione,90 giorni post ultima somministrazione
;
run;

proc print data=WORK.DESIGN_SUMMARY noobs label;
    title  "RIEPILOGO PARAMETRI DI DISEGNO DELLO STUDIO";
    var parametro valore;
    label parametro="Parametro" valore="Valore";
run;
title;


/* ── 8. RIEPILOGO IN LOG ──────────────────────────────────────────────────── */

%put NOTE: ;
%put NOTE: ══════════════════════════════════════════════════════;
%put NOTE:  RIEPILOGO CALCOLO CAMPIONE – RCT DUE GRUPPI;
%put NOTE: ══════════════════════════════════════════════════════;
%put NOTE:  Endpoint   : >=&MIN_SE SE nei &FOLLOWUP gg post ultima dose;
%put NOTE:  P controllo  (storico) : &P_CTRL;
%put NOTE:  P trattamento (storico): &P_TREAT;
%put NOTE:  P ponderata (pooled)   : &P_POOL;
%put NOTE:  Alpha (two-sided)      : &ALPHA;
%put NOTE:  Potenza                : &POWER;
%put NOTE:  Dropout                : &DROPOUT;
%put NOTE:  -> Vedi dataset WORK.SS_PRIMARY per N finale;
%put NOTE: ══════════════════════════════════════════════════════;
%put NOTE: ;


/* ── FINE PROGRAMMA ─────────────────────────────────────────────────────── */
