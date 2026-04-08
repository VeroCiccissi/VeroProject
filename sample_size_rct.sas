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
     - Gruppo controllo    (A): 44%  pazienti con ≥5 SE in 90 gg
     - Gruppo trattamento  (B): 21%  pazienti con ≥5 SE in 90 gg
     - Proporzione ponderata pooled: 30%

   MODELLO STATISTICO
   ──────────────────
   I conteggi di eventi nei due bracci sono modellati come variabili
   aleatorie con distribuzione Binomiale:

       X_ctrl  ~ Bin(n, p1)     p1 = 0.44
       X_treat ~ Bin(n, p2)     p2 = 0.21

   Stimatori di proporzione:
       p̂1 = X_ctrl  / n        E[p̂1] = p1,  Var[p̂1] = p1*(1-p1)/n
       p̂2 = X_treat / n        E[p̂2] = p2,  Var[p̂2] = p2*(1-p2)/n

   Formula N (approssimazione normale alla Binomiale – Fleiss 1981):
   ─────────────────────────────────────────────────────────────────
       n = [z_{α/2}·√(2p̄(1-p̄)) + z_β·√(p1(1-p1)+p2(1-p2))]²
           ──────────────────────────────────────────────────────
                              (p1 − p2)²

   dove p̄ = (p1+p2)/2 è la proporzione pooled sotto H0 con alloc. 1:1.
   La varianza al numeratore usa direttamente la proprietà Binomiale:
   Var[p̂] = p(1-p)/n.

   Versione con correzione per continuità (Fleiss 1981):
       n_cc = (n/4)·[1 + √(1 + 4 / (n·|p1−p2|))]²

   IPOTESI
   ────────
     H0: p1 = p2       H1: p1 ≠ p2   (two-sided)

   PARAMETRI
   ─────────
     Potenza (1−β)  : 90%     Alpha (α)  : 5%     Dropout: 10%

   Tested on SAS 9.4 / SAS Viya 4
   ============================================================================ */


/* ── 0. OPZIONI E MACRO VARIABILI ──────────────────────────────────────────── */

options nodate nonumber ls=120 ps=max;

%let P_CTRL    = 0.44;   /* X_ctrl  ~ Bin(n, 0.44)                            */
%let P_TREAT   = 0.21;   /* X_treat ~ Bin(n, 0.21)                            */
%let P_POOL    = 0.30;   /* proporzione ponderata pooled (riferimento storico) */
%let ALPHA     = 0.05;   /* livello di significatività, two-sided              */
%let POWER     = 0.90;   /* potenza desiderata                                 */
%let DROPOUT   = 0.10;   /* tasso di dropout                                   */
%let MIN_SE    = 5;      /* soglia ≥5 SE per definire l'evento                 */
%let FOLLOWUP  = 90;     /* giorni di follow-up post ultima dose               */


/* ── 1. VERIFICA PROPORZIONE PONDERATA ─────────────────────────────────────── */

data _null_;
    p1     = &P_CTRL;
    p2     = &P_TREAT;
    p_pool = &P_POOL;
    w1     = (p_pool - p2) / (p1 - p2);   /* peso implicito gruppo ctrl */
    w2     = 1 - w1;
    check  = w1*p1 + w2*p2;
    put '────────────────────────────────────────────────────────';
    put 'VERIFICA PROPORZIONE PONDERATA (studio storico)';
    put '────────────────────────────────────────────────────────';
    put '  X_ctrl  ~ Bin(n, ' p1  percent6.1 ')   peso nello storico: ' w1 percent6.1;
    put '  X_treat ~ Bin(n, ' p2  percent6.1 ')   peso nello storico: ' w2 percent6.1;
    put '  p_pool calcolata: ' check percent6.1 '   p_pool dichiarata: ' p_pool percent6.1;
    put '────────────────────────────────────────────────────────';
run;


/* ── 2. DERIVAZIONE N DALLA DISTRIBUZIONE BINOMIALE – ANALISI PRIMARIA ─────── */
/*
   Distribuzione binomiale → varianza stimatore → formula analitica per n.
   Si calcolano due versioni:
     (a) Formula Fleiss       – pooled SE sotto H0
     (b) Correzione continuità – più conservativa, raccomandata per n < 100
*/

data WORK.N_BINOMIALE_PRIMARIO;

    /* Parametri */
    p1      = &P_CTRL;     /* X_ctrl  ~ Bin(n, p1) */
    p2      = &P_TREAT;    /* X_treat ~ Bin(n, p2) */
    alpha   = &ALPHA;
    beta    = 1 - &POWER;
    dropout = &DROPOUT;

    /* Quantili della normale standard */
    z_a2    = probit(1 - alpha/2);    /* z_{α/2} = 1.9600 */
    z_b     = probit(1 - beta);       /* z_β     = 1.2816 */

    /* Proporzione pooled sotto H0 (alloc. 1:1) */
    p_bar   = (p1 + p2) / 2;

    /* Varianze binomiali dei singoli stimatori (× n) */
    v1      = p1 * (1 - p1);   /* = Var[X_ctrl]  / n  → da Bin(n,p1)  */
    v2      = p2 * (1 - p2);   /* = Var[X_treat] / n  → da Bin(n,p2)  */
    v_h0    = p_bar * (1 - p_bar);   /* varianza pooled sotto H0 */

    /* Differenza assoluta */
    delta   = abs(p1 - p2);

    /* ── (a) Formula Fleiss (1981) – approssimazione normale alla Binomiale ── */
    /*   Numeratore: usa √(2·v_h0) perché la statistica test usa la stima      */
    /*   pooled di SE sotto H0, e √(v1+v2) per la potenza sotto H1             */
    n_fleiss_raw  = ((z_a2 * sqrt(2 * v_h0) + z_b * sqrt(v1 + v2)) / delta) ** 2;
    n_fleiss      = ceil(n_fleiss_raw);

    /* ── (b) Correzione per continuità (Fleiss 1981) ── */
    /*   Aggiusta per la discretezza della Binomiale (raccomandata n < 100)    */
    n_cc_raw      = (n_fleiss / 4) * (1 + sqrt(1 + 4 / (n_fleiss * delta))) ** 2;
    n_cc          = ceil(n_cc_raw);

    /* ── Verifica potenza effettiva a n_fleiss (approx. normale) ── */
    se_h1         = sqrt((v1 + v2) / n_fleiss);
    z_eff         = delta / se_h1 - z_a2;
    power_eff     = probnorm(z_eff);

    /* ── Correzione dropout ── */
    n_fleiss_adj  = ceil(n_fleiss / (1 - dropout));
    n_totale_fl   = n_fleiss_adj  * 2;
    n_cc_adj      = ceil(n_cc     / (1 - dropout));
    n_totale_cc   = n_cc_adj      * 2;

    label
        p1            = "P Controllo  [X_ctrl  ~ Bin(n,p1)]"
        p2            = "P Trattamento [X_treat ~ Bin(n,p2)]"
        p_bar         = "P Pooled H0  [(p1+p2)/2]"
        v1            = "Var Binomiale ctrl  = p1(1-p1)"
        v2            = "Var Binomiale treat = p2(1-p2)"
        delta         = "|p1 - p2|"
        z_a2          = "z_{α/2}"
        z_b           = "z_{β}"
        n_fleiss      = "N/gruppo – Fleiss"
        n_cc          = "N/gruppo – Fleiss + cont. correction"
        power_eff     = "Potenza effettiva a n_fleiss (verifica)"
        n_fleiss_adj  = "N/gruppo Fleiss + dropout 10%"
        n_totale_fl   = "N Totale Fleiss + dropout 10%"
        n_cc_adj      = "N/gruppo CC + dropout 10%"
        n_totale_cc   = "N Totale CC + dropout 10%";
run;

proc print data=WORK.N_BINOMIALE_PRIMARIO noobs label;
    title  "ANALISI PRIMARIA – N da Distribuzione Binomiale";
    title2 "Endpoint: ≥&MIN_SE SE in &FOLLOWUP gg | X_ctrl~Bin(n,0.44) vs X_treat~Bin(n,0.21)";

    var p1 p2 p_bar delta z_a2 z_b v1 v2
        n_fleiss power_eff n_fleiss_adj n_totale_fl
        n_cc     n_cc_adj  n_totale_cc;

    format p1 p2 p_bar          percent6.1
           delta                 6.4
           z_a2 z_b              6.4
           v1 v2                 7.4
           n_fleiss n_cc         6.0
           power_eff             percent6.1
           n_fleiss_adj n_cc_adj
           n_totale_fl n_totale_cc  6.0;
run;
title;


/* ── 3. VERIFICA CON TEST ESATTO (Fisher / Binomiale esatta) ────────────────── */
/*
   PROC POWER TEST=FISH usa la distribuzione Binomiale esatta (non
   l'approssimazione normale) per calcolare la potenza e il campione.
   Il confronto con la formula di Fleiss valida la stima analitica.
*/

proc power;
    title  "VERIFICA – Test Esatto Binomiale (Fisher)";
    title2 "Confronto con formula Fleiss | X_ctrl~Bin(n,0.44) vs X_treat~Bin(n,0.21)";

    twosamplefreq test = fish              /* distribuzione Binomiale esatta   */
        groupproportions = (&P_CTRL &P_TREAT)
        alpha            = &ALPHA
        power            = &POWER
        sides            = 2
        npergroup        = .;

    ods output output = _power_exact;
run;
title;

data WORK.N_ESATTO;
    set _power_exact;
    n_adj    = ceil(NPerGroup / (1 - &DROPOUT));
    n_totale = n_adj * 2;
    fonte    = "PROC POWER TEST=FISH (binomiale esatta)";
    label
        NPerGroup = "N/gruppo (esatto)"
        n_adj     = "N/gruppo + dropout 10%"
        n_totale  = "N Totale + dropout 10%";
    keep fonte NPerGroup NominalPower n_adj n_totale;
run;

proc print data=WORK.N_ESATTO noobs label;
    title  "VERIFICA Test Binomiale Esatto + Dropout 10%";
    format NominalPower percent6.1
           NPerGroup n_adj n_totale 6.0;
run;
title;


/* ── 4. RIEPILOGO COMPARATIVO: FORMULA vs TEST ESATTO ──────────────────────── */

data WORK.CONFRONTO;
    length metodo $55;

    set WORK.N_BINOMIALE_PRIMARIO (
            keep=n_fleiss n_fleiss_adj n_totale_fl power_eff
            rename=(n_fleiss=n_pgruppo n_fleiss_adj=n_adj n_totale_fl=n_totale
                    power_eff=potenza_effettiva))
        WORK.N_ESATTO (
            rename=(NPerGroup=n_pgruppo NominalPower=potenza_effettiva));

    if _n_ = 1 then metodo = "Formula Fleiss (normal approx binomiale)";
    else            metodo = "PROC POWER TEST=FISH (binomiale esatta)";

    label
        metodo            = "Metodo"
        n_pgruppo         = "N/gruppo (no dropout)"
        n_adj             = "N/gruppo (dropout +10%)"
        n_totale          = "N Totale (dropout +10%)"
        potenza_effettiva = "Potenza";
run;

proc print data=WORK.CONFRONTO noobs label;
    title  "CONFRONTO METODI – Formula Fleiss vs Binomiale Esatta";
    title2 "Endpoint: ≥&MIN_SE SE in &FOLLOWUP gg | Dropout 10% | p_ctrl=44% p_treat=21%";
    var metodo n_pgruppo n_adj n_totale potenza_effettiva;
    format potenza_effettiva percent6.1
           n_pgruppo n_adj n_totale 6.0;
run;
title;


/* ── 5. ANALISI DI SENSIBILITÀ – FORMULA BINOMIALE ─────────────────────────── */
/*
   Si mantiene il gruppo controllo a p1=44% [X_ctrl ~ Bin(n, 0.44)].
   Si varia p2 (proporzione trattamento) da 0.15 a 0.35.
   N calcolato con formula Fleiss su varianza binomiale.
   Il valore storico p2=21% è marcato con (*).
*/

data WORK.SENSIBILITA_BINOMIALE;

    p1      = &P_CTRL;
    alpha   = &ALPHA;
    z_a2    = probit(1 - alpha/2);
    z_b     = probit(&POWER);
    dropout = &DROPOUT;

    do p2 = 0.15, 0.21, 0.25, 0.30, 0.35;

        /* Varianze binomiali */
        v1    = p1 * (1 - p1);
        v2    = p2 * (1 - p2);
        p_bar = (p1 + p2) / 2;
        delta = abs(p1 - p2);

        /* Formula Fleiss */
        n_raw  = ((z_a2 * sqrt(2 * p_bar*(1-p_bar))
                 + z_b  * sqrt(v1 + v2)) / delta) ** 2;
        n_base = ceil(n_raw);

        /* Correzione continuità */
        n_cc   = ceil((n_base/4) * (1 + sqrt(1 + 4/(n_base*delta)))**2);

        /* Potenza effettiva (verifica) */
        se_h1       = sqrt((v1+v2) / n_base);
        power_eff   = probnorm(delta/se_h1 - z_a2);

        /* Dropout */
        n_adj   = ceil(n_base / (1 - dropout));
        n_total = n_adj * 2;

        /* Flag valore storico */
        storico = (abs(p2 - &P_TREAT) < 0.001);

        output;
    end;

    label
        p1          = "P Ctrl  [Bin(n,p1)]"
        p2          = "P Treat [Bin(n,p2)]"
        delta       = "|p1-p2| (pp)"
        n_base      = "N/grp Fleiss"
        n_cc        = "N/grp CC"
        power_eff   = "Potenza effettiva"
        n_adj       = "N/grp + dropout"
        n_total     = "N Totale + dropout"
        storico     = "(*) Valore storico";

    keep p1 p2 delta n_base n_cc power_eff n_adj n_total storico;
run;

proc print data=WORK.SENSIBILITA_BINOMIALE noobs label;
    title  "SENSIBILITÀ – Variazione P Trattamento | Formula Binomiale";
    title2 "X_ctrl~Bin(n,0.44) fisso | X_treat~Bin(n,p2) variabile | Dropout 10%";
    title3 "(*) = scenario storico principale (p_treat=21%)";

    var p1 p2 delta n_base n_cc power_eff n_adj n_total storico;

    format p1 p2           percent6.1
           delta            5.2
           n_base n_cc      6.0
           power_eff        percent6.1
           n_adj n_total    6.0
           storico          1.0;
run;
title;


/* ── 6. CURVA DI POTENZA – DISTRIBUZIONE BINOMIALE ESATTA ──────────────────── */

proc power plotonly;
    title  "Curva di Potenza – Test Binomiale Esatto (Fisher)";
    title2 "X_ctrl~Bin(n,0.44) | Variazione p_treat | Alpha=5% | Two-sided";

    twosamplefreq test = fish
        refproportion  = &P_CTRL
        proportiondiff = -0.09 -0.14 -0.19 -0.23 -0.29
        /* p_treat:        0.35   0.30   0.25   0.21*  0.15 */
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


/* ── 7. SCENARIO ALTERNATIVO: PROPORZIONE UNICA PONDERATA p=0.30 ────────────── */
/*
   Anziché confrontare due gruppi (44% vs 21%), si usa la singola proporzione
   ponderata p=0.30 derivata dalla media pesata dei due bracci storici.

   MODELLO ONE-SAMPLE:
       X ~ Bin(n, p)     con p = 0.30   (proporzione pooled ponderata)

   TEST:
       H0: p = p0   (valore null di riferimento)
       H1: p = 0.30 (proporzione pooled osservata ≠ p0)

   FORMULA (approssimazione normale alla Binomiale):
       n = [z_{α/2}·√(p0(1-p0))  +  z_β·√(p1(1-p1))]²
           ──────────────────────────────────────────────
                          (p1 − p0)²

   dove p1=0.30 (H1) e p0 è il valore null da specificare.

   Valori null considerati:
     p0=0.44  → testare che p_pool < p_ctrl  (efficacia rispetto al controllo)
     p0=0.40  → scenario conservativo
     p0=0.21  → testare che p_pool > p_treat (sicurezza rispetto al trattamento)
     p0=0.50  → testare che l'incidenza sia sotto la soglia del 50%
*/

/* ── 7a. Formula one-sample binomiale – DATA step ── */

data WORK.N_ONESAMPLE_BINOM;

    p1      = &P_POOL;       /* X ~ Bin(n, 0.30) sotto H1                  */
    alpha   = &ALPHA;
    dropout = &DROPOUT;

    z_a2    = probit(1 - alpha/2);    /* z_{α/2} = 1.9600 */
    z_b     = probit(&POWER);         /* z_{β}   = 1.2816 */
    v1      = p1 * (1 - p1);          /* Var binomiale sotto H1: 0.30*0.70  */

    do p0 = 0.44, 0.40, 0.21, 0.50;

        delta = abs(p1 - p0);
        v0    = p0 * (1 - p0);        /* Var binomiale sotto H0: p0*(1-p0)  */

        /* Formula one-sample (approssimazione normale alla Binomiale) */
        n_raw  = ((z_a2 * sqrt(v0) + z_b * sqrt(v1)) / delta) ** 2;
        n_base = ceil(n_raw);

        /* Correzione continuità */
        n_cc   = ceil(n_base + 1 / (2 * n_base * delta**2));

        /* Potenza effettiva a n_base */
        se_h1     = sqrt(v1 / n_base);
        power_eff = probnorm(delta / se_h1 - z_a2);

        /* Dropout */
        n_adj    = ceil(n_base / (1 - dropout));

        /* Etichetta interpretativa */
        if      p0 = 0.44 then interp = "p=0.30 < p_ctrl(0.44): test efficacia";
        else if p0 = 0.40 then interp = "p=0.30 < 0.40: scenario conservativo";
        else if p0 = 0.21 then interp = "p=0.30 > p_treat(0.21): test sicurezza";
        else if p0 = 0.50 then interp = "p=0.30 < 0.50: incidenza sotto soglia";

        output;
    end;

    label
        p1        = "P Osservata H1  [X~Bin(n,0.30)]"
        p0        = "P Null H0"
        delta     = "|p1 - p0|"
        v0        = "Var Bin H0 = p0(1-p0)"
        v1        = "Var Bin H1 = p1(1-p1)"
        n_base    = "n (formula one-sample)"
        n_cc      = "n (cont. correction)"
        power_eff = "Potenza effettiva"
        n_adj     = "n + dropout 10%"
        interp    = "Interpretazione clinica";

    keep p1 p0 delta v0 v1 n_base n_cc power_eff n_adj interp;
run;

proc print data=WORK.N_ONESAMPLE_BINOM noobs label;
    title  "SCENARIO ALTERNATIVO – Test One-Sample, X~Bin(n, 0.30)";
    title2 "Proporzione unica ponderata p=0.30 | Power=90% | Alpha=5% | Dropout=10%";
    title3 "Endpoint: ≥&MIN_SE SE in &FOLLOWUP gg post ultima dose";
    var p1 p0 delta v0 v1 n_base n_cc power_eff n_adj interp;
    format p1 p0       percent6.1
           delta v0 v1 7.4
           n_base n_cc 6.0
           power_eff   percent6.1
           n_adj       6.0;
run;
title;


/* ── 7b. Verifica con PROC POWER ONESAMPLEFREQ (binomiale esatta) ── */

proc power;
    title  "One-Sample – Test Binomiale Esatto su p=0.30";
    title2 "Verifica formula | Alpha=5% | Power=90% | Two-sided";

    onesamplefreq test = z
        proportion     = &P_POOL              /* p1 = 0.30 sotto H1           */
        nullproportion = 0.44 0.40 0.21 0.50  /* diversi valori null          */
        alpha          = &ALPHA
        power          = &POWER
        sides          = 2
        ntotal         = .;

    ods output output = _power_1s;
run;
title;

data WORK.N_ONESAMPLE_ADJ;
    set _power_1s;
    n_adj    = ceil(NTotal / (1 - &DROPOUT));
    label
        NullProportion = "P Null H0"
        Proportion     = "P H1 (pooled)"
        NominalPower   = "Potenza"
        NTotal         = "n (PROC POWER)"
        n_adj          = "n + dropout 10%";
    keep NullProportion Proportion NominalPower NTotal n_adj;
run;

proc print data=WORK.N_ONESAMPLE_ADJ noobs label;
    title  "One-Sample – Verifica PROC POWER + Dropout 10%";
    format NullProportion Proportion percent6.1
           NominalPower percent6.1
           NTotal n_adj 6.0;
run;
title;


/* ── 7c. Confronto one-sample vs two-sample ── */
/*
   Riepilogo side-by-side dei due approcci per lo stesso studio.
   Il confronto guida la scelta del disegno più appropriato.
*/

data WORK.CONFRONTO_DISEGNI;
    length disegno $55 ipotesi $55;

    /* Due campioni (primario) – da WORK.N_BINOMIALE_PRIMARIO */
    disegno  = "Two-sample (p_ctrl=44% vs p_treat=21%)";
    ipotesi  = "H0: p_ctrl=p_treat  H1: p_ctrl≠p_treat";
    n_pgruppo = .;
    n_totale  = .;

    /* Leggi i valori dal dataset primario */
    set WORK.N_BINOMIALE_PRIMARIO
            (keep = n_fleiss n_fleiss_adj n_totale_fl
             rename=(n_fleiss=n_pgruppo n_fleiss_adj=_adj n_totale_fl=n_totale));
    output; stop;
run;

proc sql noprint;
    /* Aggiunge la riga one-sample per p0=0.44 (confronto più rilevante) */
    insert into WORK.CONFRONTO_DISEGNI
    select "One-sample (p_pool=30% vs H0: p=44%)" as disegno,
           "H0: p=0.44  H1: p=0.30" as ipotesi,
           . as n_pgruppo,
           n_adj as n_totale
    from WORK.N_ONESAMPLE_ADJ
    where NullProportion = 0.44;
quit;

proc print data=WORK.CONFRONTO_DISEGNI noobs label;
    title  "CONFRONTO DISEGNI – Two-sample vs One-sample";
    title2 "Stessa potenza 90%, alpha 5%, dropout 10%";
    var disegno ipotesi n_pgruppo n_totale;
    label disegno   = "Disegno"
          ipotesi   = "Ipotesi statistica"
          n_pgruppo = "N per gruppo"
          n_totale  = "N totale";
    format n_pgruppo n_totale 6.0;
run;
title;


/* ── 8. RIEPILOGO DESIGN IN LOG ─────────────────────────────────────────────── */

%put NOTE: ;
%put NOTE: ══════════════════════════════════════════════════════════════════;
%put NOTE:  SAMPLE SIZE – DISTRIBUZIONE BINOMIALE;
%put NOTE: ══════════════════════════════════════════════════════════════════;
%put NOTE:  Endpoint   : >= &MIN_SE SE nei &FOLLOWUP gg post ultima dose;
%put NOTE: ;
%put NOTE:  [DISEGNO A – DUE CAMPIONI] (analisi primaria);
%put NOTE:   X_ctrl~Bin(n,&P_CTRL)  vs  X_treat~Bin(n,&P_TREAT);
%put NOTE:   Metodo: Formula Fleiss + PROC POWER TEST=FISH;
%put NOTE:   -> WORK.N_BINOMIALE_PRIMARIO | WORK.CONFRONTO;
%put NOTE: ;
%put NOTE:  [DISEGNO B – ONE-SAMPLE] (scenario alternativo);
%put NOTE:   X~Bin(n,&P_POOL)  |  p_pool=&P_POOL da media ponderata;
%put NOTE:   H0: p=p0 (valori testati: 0.44 / 0.40 / 0.21 / 0.50);
%put NOTE:   Metodo: Formula one-sample + PROC POWER ONESAMPLEFREQ;
%put NOTE:   -> WORK.N_ONESAMPLE_BINOM | WORK.N_ONESAMPLE_ADJ;
%put NOTE: ;
%put NOTE:  Parametri comuni: alpha=&ALPHA | power=&POWER | dropout=&DROPOUT;
%put NOTE: ══════════════════════════════════════════════════════════════════;
%put NOTE: ;


/* ── FINE PROGRAMMA ─────────────────────────────────────────────────────── */
