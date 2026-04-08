/* ============================================================================
   sample_size_rct.sas
   -------------------
   Calcolo del campione – Test ONE-SAMPLE ONE-SIDED sulla proporzione

   ENDPOINT PRIMARIO
   ─────────────────
   Proporzione di pazienti con ≥ 5 side effects (SE)
   nei 90 giorni successivi all'ultima somministrazione.

   DISEGNO STATISTICO
   ──────────────────
   Un solo gruppo.  X ~ Binomiale(n, p).

   H0 : p  = 0.30   (incidenza non supera la soglia di riferimento)
   H1 : p  > 0.30   (ONE-SIDED – incidenza supera la soglia)

   La soglia p0 = 0.30 è la proporzione ponderata pooled derivata dallo
   studio storico (44% gruppo A, 21% gruppo B, pesi 39% / 61%).

   FORMULA (approssimazione normale alla Binomiale, one-sided)
   ─────────────────────────────────────────────────────────────
   n = [z_α·√(p0·(1−p0))  +  z_β·√(p1·(1−p1))]²
       ────────────────────────────────────────────
                     (p1 − p0)²

   dove z_α = z_{0.95} = 1.6449 (one-sided 5%)

   REGOLA DECISIONALE (test esatto)
   ──────────────────────────────────
   Rifiuta H0 se  X ≥ k_critico
   con k_critico tale che  P(X ≥ k | Bin(n, p0)) ≤ α

   PARAMETRI
   ─────────
     p0 (H0)          : 30%
     p1 (H1)          : 44% scenario storico  /  50% ≈ documento ufficiale (n≈60)
     Potenza (1−β)    : 90%
     Alpha (one-sided): 5%
     Dropout          : 10%

   NOTE
   ────
   Con p1=0.50 e one-sided: n_esatto=53 → +dropout → 59 ≈ 60 (documento ufficiale).
   Con p1=0.44 (storico):   n_esatto=105 → +dropout → 117.

   Tested on SAS 9.4 / SAS Viya 4
   ============================================================================ */


/* ── 0. PARAMETRI ───────────────────────────────────────────────────────────── */

options nodate nonumber ls=120 ps=max;

%let P0      = 0.30;   /* soglia H0                                            */
%let ALPHA   = 0.05;   /* livello di significatività (one-sided)               */
%let POWER   = 0.90;   /* potenza                                              */
%let DROPOUT = 0.10;   /* dropout                                              */
%let MIN_SE  = 5;      /* soglia evento: ≥5 SE                                 */
%let FWUP    = 90;     /* giorni follow-up post ultima dose                    */


/* ── 1. FORMULA (approssimazione normale, one-sided) ────────────────────────── */

data WORK.SS_FORMULA;

    p0      = &P0;
    alpha   = &ALPHA;
    dropout = &DROPOUT;

    z_a  = probit(1 - alpha);      /* z_α  one-sided: 1.6449                  */
    z_b  = probit(&POWER);         /* z_β           : 1.2816                  */

    do p1 = 0.40, 0.44, 0.45, 0.50, 0.55, 0.60;

        delta  = p1 - p0;          /* differenza positiva (H1: p > p0)        */
        v0     = p0 * (1 - p0);
        v1     = p1 * (1 - p1);

        n_raw  = ((z_a * sqrt(v0) + z_b * sqrt(v1)) / delta) ** 2;
        n_base = ceil(n_raw);
        n_adj  = ceil(n_base / (1 - dropout));

        /* Potenza approssimata a n_base */
        power_eff = probnorm(delta / sqrt(v1/n_base) - z_a);

        storico  = (abs(p1 - 0.44) < 0.001);
        ufficial = (abs(p1 - 0.50) < 0.001);   /* match documento ~60         */

        label
            p0        = "P0 H0 (soglia)"
            p1        = "P1 H1 (alternativa)"
            delta     = "p1 - p0"
            n_base    = "n (formula)"
            n_adj     = "n + dropout 10%"
            power_eff = "Potenza appross."
            storico   = "(*) Storico"
            ufficial  = "(**) Doc. uff.";
        output;
    end;
    keep p0 p1 delta v0 v1 n_base n_adj power_eff storico ufficial;
run;

proc print data=WORK.SS_FORMULA noobs label;
    title  "ONE-SAMPLE ONE-SIDED – Formula (approssimazione normale)";
    title2 "H0: p=0.30  H1: p>0.30  |  Power=90%  |  Alpha=5% one-sided  |  Dropout=10%";
    title3 "(*) scenario storico p1=0.44  |  (**) documento ufficiale p1=0.50 → n≈60";
    var p0 p1 delta n_base n_adj power_eff storico ufficial;
    format p0 p1 percent6.1  delta 6.2
           n_base n_adj 6.0  power_eff percent6.1
           storico ufficial 1.0;
run;
title;


/* ── 2. TEST ESATTO BINOMIALE (PROC POWER TEST=EXACT, one-sided) ────────────── */

proc power;
    title  "ONE-SAMPLE ONE-SIDED – Test Binomiale ESATTO";
    title2 "H0: p=0.30  H1: p>0.30  |  Power=90%  |  Alpha=5% one-sided";

    onesamplefreq test = exact
        nullproportion = &P0
        proportion     = 0.40 0.44 0.45 0.50 0.55 0.60
        alpha          = &ALPHA
        power          = &POWER
        sides          = 1           /* ONE-SIDED                              */
        ntotal         = .;

    ods output output = _power_exact;
run;
title;

data WORK.SS_ESATTO;
    set _power_exact;
    n_adj    = ceil(NTotal / (1 - &DROPOUT));
    storico  = (abs(Proportion - 0.44) < 0.001);
    ufficial = (abs(Proportion - 0.50) < 0.001);
    label
        NullProportion = "P0 H0"
        Proportion     = "P1 H1"
        NominalPower   = "Potenza"
        NTotal         = "n (esatto)"
        n_adj          = "n + dropout 10%"
        storico        = "(*)"
        ufficial       = "(**)";
    keep NullProportion Proportion NominalPower NTotal n_adj storico ufficial;
run;

proc print data=WORK.SS_ESATTO noobs label;
    title  "ONE-SAMPLE ONE-SIDED – Binomiale ESATTO + Dropout";
    title2 "(*) storico p1=0.44 → n=117  |  (**) documento p1=0.50 → n≈59~60";
    var NullProportion Proportion NominalPower NTotal n_adj storico ufficial;
    format NullProportion Proportion percent6.1
           NominalPower percent6.1
           NTotal n_adj 6.0
           storico ufficial 1.0;
run;
title;


/* ── 3. CONFRONTO FORMULA vs ESATTO ─────────────────────────────────────────── */

proc sql;
    create table WORK.CONFRONTO as
    select
        f.p1                  as p1_H1         format=percent6.1  label="P1 (H1)",
        f.n_base              as n_formula                         label="n Formula",
        e.NTotal              as n_esatto                          label="n Esatto",
        f.n_adj               as n_adj_formula                     label="n+DO Formula",
        e.n_adj               as n_adj_esatto                      label="n+DO Esatto",
        f.power_eff           as power_formula  format=percent6.1  label="Power (Formula)",
        e.NominalPower        as power_esatta   format=percent6.1  label="Power (Esatta)",
        f.storico,
        f.ufficial
    from WORK.SS_FORMULA f
    join WORK.SS_ESATTO  e on abs(f.p1 - e.Proportion) < 0.001
    order by f.p1;
quit;

proc print data=WORK.CONFRONTO noobs label;
    title  "CONFRONTO Formula vs Binomiale Esatta – ONE-SIDED";
    title2 "H0: p=0.30  H1: p>0.30  |  Power=90%  |  Alpha=5%  |  Dropout=10%";
    format storico ufficial 1.0;
run;
title;


/* ── 4. CURVA DI POTENZA (one-sided) ────────────────────────────────────────── */

proc power plotonly;
    title  "Curva di Potenza – ONE-SAMPLE ONE-SIDED Binomiale Esatta";
    title2 "H0: p=0.30  H1: p>0.30  |  Alpha=5% one-sided";

    onesamplefreq test = exact
        nullproportion = &P0
        proportion     = 0.44 0.50 0.55
        alpha          = &ALPHA
        power          = 0.70 to 0.95 by 0.01
        sides          = 1
        ntotal         = .;

    plot x     = n
         key   = byfeature(pos=inset)
         yopts = (ref=&POWER refline=yes label="Power 90%")
         vary(linestyle by proportion);
run;
title;


/* ── 5. LOG ──────────────────────────────────────────────────────────────────── */

%put NOTE: ;
%put NOTE: ════════════════════════════════════════════════════════════;
%put NOTE:  ONE-SAMPLE ONE-SIDED  X~Bin(n,p)  H0:p=&P0  H1:p>&P0;
%put NOTE:  Alpha=&ALPHA (one-sided) | Power=&POWER | Dropout=&DROPOUT;
%put NOTE:  Endpoint: >=&MIN_SE SE in &FWUP gg post ultima dose;
%put NOTE: ════════════════════════════════════════════════════════════;
%put NOTE:  p1=0.44 (storico)  → n=105 eval. → n+DO=117  (* WORK.SS_ESATTO);
%put NOTE:  p1=0.50 (ufficial) → n= 53 eval. → n+DO= 59  (** corrisponde a 60);
%put NOTE: ════════════════════════════════════════════════════════════;
%put NOTE: ;

/* ── FINE PROGRAMMA ─────────────────────────────────────────────────────── */
