/* ============================================================================
   sample_size_rct.sas
   -------------------
   Calcolo del campione – Test ONE-SAMPLE a due code sulla proporzione

   ENDPOINT PRIMARIO
   ─────────────────
   Proporzione di pazienti con ≥ 5 side effects (SE)
   nei 90 giorni successivi all'ultima somministrazione.

   DISEGNO STATISTICO
   ──────────────────
   Un solo gruppo.  X ~ Binomiale(n, p).

   H0 : p  = 0.30   (incidenza non diversa dalla soglia di riferimento)
   H1 : p ≠ 0.30   (two-sided)

   La soglia p0 = 0.30 è la proporzione ponderata pooled derivata dallo
   studio storico (44% gruppo A, 21% gruppo B, pesi 39% / 61%).

   FORMULA (approssimazione normale alla Binomiale)
   ─────────────────────────────────────────────────
   n = [z_{α/2}·√(p0·(1−p0))  +  z_β·√(p1·(1−p1))]²
       ────────────────────────────────────────────────
                        (p1 − p0)²

   PARAMETRI
   ─────────
     p0 (H0 – soglia)        : 30%
     p1 (H1 – alternativa)   : 44% scenario storico  /  50% scenario ufficiale
     Potenza (1−β)           : 90%
     Alpha (two-sided)       : 5%
     Dropout                 : 10%

   Tested on SAS 9.4 / SAS Viya 4
   ============================================================================ */


/* ── 0. PARAMETRI ───────────────────────────────────────────────────────────── */

options nodate nonumber ls=120 ps=max;

%let P0      = 0.30;   /* soglia H0  – proporzione di riferimento             */
%let ALPHA   = 0.05;   /* livello di significatività (two-sided)              */
%let POWER   = 0.90;   /* potenza                                             */
%let DROPOUT = 0.10;   /* tasso di dropout                                    */
%let MIN_SE  = 5;      /* soglia evento: ≥5 SE                                */
%let FWUP    = 90;     /* giorni di follow-up post ultima dose                */


/* ── 1. CALCOLO N CON FORMULA (approx. normale) ─────────────────────────────── */
/*
   n = [z_{α/2}·√(p0(1-p0))  +  z_β·√(p1(1-p1))]²  /  (p1-p0)²
   Varianza binomiale: Var[p̂] = p(1-p)/n
*/

data WORK.SS_FORMULA;

    p0      = &P0;
    alpha   = &ALPHA;
    dropout = &DROPOUT;

    z_a2    = probit(1 - alpha/2);   /* z_{α/2} = 1.9600 */
    z_b     = probit(&POWER);        /* z_{β}   = 1.2816 */

    do p1 = 0.40, 0.44, 0.45, 0.50, 0.55, 0.60;

        delta = abs(p1 - p0);
        v0    = p0 * (1 - p0);      /* Var binomiale sotto H0 */
        v1    = p1 * (1 - p1);      /* Var binomiale sotto H1 */

        n_raw  = ((z_a2 * sqrt(v0) + z_b * sqrt(v1)) / delta) ** 2;
        n_base = ceil(n_raw);
        n_adj  = ceil(n_base / (1 - dropout));

        /* Potenza effettiva a n_base */
        power_eff = probnorm(delta / sqrt(v1/n_base) - z_a2)
                  + probnorm(-delta / sqrt(v1/n_base) - z_a2);

        storico  = (abs(p1 - 0.44) < 0.001);   /* scenario storico */
        ufficial = (abs(p1 - 0.50) < 0.001);   /* ~n=60 ufficiale  */

        label
            p0        = "P0 (H0)"
            p1        = "P1 (H1)"
            delta     = "|p1-p0| (pp)"
            v0        = "Var Bin H0"
            v1        = "Var Bin H1"
            n_base    = "n (formula)"
            n_adj     = "n + dropout 10%"
            power_eff = "Potenza effettiva"
            storico   = "(*) Storico"
            ufficial  = "(**) Ufficiale";

        output;
    end;

    keep p0 p1 delta v0 v1 n_base n_adj power_eff storico ufficial;
run;

proc print data=WORK.SS_FORMULA noobs label;
    title  "ONE-SAMPLE – Formula (approssimazione normale alla Binomiale)";
    title2 "H0: p=0.30  H1: p≠0.30  |  Power=90%  |  Alpha=5% two-sided  |  Dropout=10%";
    title3 "Endpoint: ≥&MIN_SE SE in &FWUP gg post ultima dose";
    var p0 p1 delta v0 v1 n_base n_adj power_eff storico ufficial;
    format p0 p1      percent6.1
           delta v0 v1 7.4
           n_base n_adj 6.0
           power_eff  percent6.1
           storico ufficial 1.0;
run;
title;


/* ── 2. CALCOLO N CON TEST ESATTO BINOMIALE ─────────────────────────────────── */
/*
   PROC POWER TEST=EXACT usa la distribuzione Binomiale esatta
   (non l'approssimazione normale) per calcolare potenza e n.
   Regola decisionale: rifiuta H0 se X ≤ k_low OPPURE X ≥ k_high.
*/

proc power;
    title  "ONE-SAMPLE – Test Binomiale ESATTO (PROC POWER TEST=EXACT)";
    title2 "H0: p=0.30  H1: p≠0.30  |  Power=90%  |  Alpha=5% two-sided";

    onesamplefreq test = exact
        nullproportion = &P0
        proportion     = 0.40 0.44 0.45 0.50 0.55 0.60
        alpha          = &ALPHA
        power          = &POWER
        sides          = 2
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
        NullProportion = "P0 (H0)"
        Proportion     = "P1 (H1)"
        NominalPower   = "Potenza"
        NTotal         = "n (esatto)"
        n_adj          = "n + dropout 10%"
        storico        = "(*) Storico"
        ufficial       = "(**) Ufficiale";
    keep NullProportion Proportion NominalPower NTotal n_adj storico ufficial;
run;

proc print data=WORK.SS_ESATTO noobs label;
    title  "ONE-SAMPLE – Test Binomiale ESATTO + Dropout 10%";
    title2 "(*) scenario storico p1=0.44  |  (**) ~n=60 documento ufficiale p1=0.50";
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
        f.p1                        as P1_H1     format=percent6.1,
        f.n_base                    as n_formula,
        e.NTotal                    as n_esatto,
        f.n_adj                     as n_adj_formula,
        e.n_adj                     as n_adj_esatto,
        f.power_eff                 as power_formula format=percent6.1,
        f.storico,
        f.ufficial
    from WORK.SS_FORMULA   as f
    join WORK.SS_ESATTO    as e
      on abs(f.p1 - e.Proportion) < 0.001
    order by f.p1;
quit;

proc print data=WORK.CONFRONTO noobs label;
    title  "CONFRONTO Formula vs Esatto";
    title2 "H0: p=0.30  H1: p≠0.30  |  Power=90%  |  Alpha=5%  |  Dropout=10%";
    label P1_H1         = "P1 (H1)"
          n_formula     = "n Formula"
          n_esatto      = "n Esatto"
          n_adj_formula = "n+DO Formula"
          n_adj_esatto  = "n+DO Esatto"
          power_formula = "Power (Formula)"
          storico       = "(*)"
          ufficial      = "(**)";
    format P1_H1 percent6.1 power_formula percent6.1 storico ufficial 1.0;
run;
title;


/* ── 4. CURVA DI POTENZA ─────────────────────────────────────────────────────── */

proc power plotonly;
    title  "Curva di Potenza – ONE-SAMPLE Binomiale Esatta";
    title2 "H0: p=0.30  |  Alpha=5% two-sided  |  X ~ Bin(n, p)";

    onesamplefreq test = exact
        nullproportion = &P0
        proportion     = 0.44 0.50 0.55
        alpha          = &ALPHA
        power          = 0.70 to 0.95 by 0.01
        sides          = 2
        ntotal         = .;

    plot x     = n
         key   = byfeature(pos=inset)
         yopts = (ref=&POWER refline=yes label="Power 90%")
         vary(linestyle by proportion);
run;
title;


/* ── 5. RIEPILOGO IN LOG ─────────────────────────────────────────────────────── */

%put NOTE: ;
%put NOTE: ══════════════════════════════════════════════════════════════════;
%put NOTE:  ONE-SAMPLE BINOMIALE – TEST A DUE CODE;
%put NOTE:  X ~ Bin(n, p)   H0: p=&P0   H1: p≠&P0;
%put NOTE:  Endpoint: ≥&MIN_SE SE in &FWUP gg | Alpha=&ALPHA | Power=&POWER | DO=&DROPOUT;
%put NOTE: ══════════════════════════════════════════════════════════════════;
%put NOTE:  Risultati (n binomiale ESATTO + dropout):;
%put NOTE:  p1=0.44 (storico)  -> vedere WORK.SS_ESATTO;
%put NOTE:  p1=0.50 (ufficial) -> vedere WORK.SS_ESATTO;
%put NOTE:  Confronto formula vs esatto -> WORK.CONFRONTO;
%put NOTE: ══════════════════════════════════════════════════════════════════;
%put NOTE: ;


/* ── FINE PROGRAMMA ─────────────────────────────────────────────────────── */
