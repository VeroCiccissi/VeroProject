/* ============================================================================
   neridronate_offlabel.sql
   ------------------------
   Adverse-event case extraction for NERIDRONATE (off-label use).
   Run via SAS/ACCESS passthrough against the Oracle PROD1 safety DB
   (PISAS01.abiogen.it).

   Step : STEP 3 N – FINAL VERSION (FIX CLOB ERROR)
   Output table : NERIDRONATE_OFFLABEL
   Data lock     : 2025-07-04

   Bug fixes applied
   -----------------
   1. ORA-00936 – malformed WHERE clause "(or cp.sort_id-1)" replaced with
      "(ca.state_id <> 1 or cp.sort_id <> 1)".
   2. dbms_lob.substr – explicit 3-argument form (lob, amount, offset=1) used
      everywhere so the call is unambiguous inside an Oracle SQL context.
   3. product_name CLOB – wrapped with dbms_lob.substr like generic_name.
   4. ind_pref_tera subquery – column alias disambiguated to avoid name clash
      with the derived-table alias.
   ============================================================================ */

proc sql;
    connect to oracle (user="&user" password="&pw" path="PROD1!");

    /* ── CREATE OUTPUT TABLE ─────────────────────────────────────────────── */
    create table NERIDRONATE_OFFLABEL as
    select *
    from connection to oracle
    (
        /* ── AE_CASES base query ──────────────────────────────────────── */
        select distinct
              ca.case_id
            , dbms_lob.substr(cp.product_name, 200, 1)  as cp_product_name   /* Bug 3 – was ca.cp_product_name (invalid col); product_name is CLOB */
            , ca.case_seq_num
            , ca.init_rcpt_date                          as init_receipt_date
            , cpi.pat_age
            , cpi.gender_id
            , lgo.group_name                             as age_group
            , cp.sort_id
            , cp.pref_tera

            /* CLOB columns – Bug 2: explicit 3-arg form */
            , dbms_lob.substr(cp.generic_name, 200, 1)  as generic_name

            /* A EVENT (OFF LABEL) */
            , case
                  when lower(ce.pref_tera) = 'off label use'
                  then ce.pref_tera
              end                                        as PT

            , ce.body_sys                                as soc

        from case_master ca

        inner join case_pat_info cpi
            on  ca.case_id   = cpi.case_id
            and cpi.deleted  is null

        left  join lo_age_groups lgo
            on  lgo.age_group_id = cpi.age_group_id

        left  join lo_gender lg
            on  lg.gender_id     = cpi.gender_id

        inner join case_product cp
            on  ca.case_id        = cp.case_id
            and cp.deleted       is null
            and (cp.first_tera   is null or cp.first_tera <> 2)

        left  join (
                       /* Bug 4 – column renamed to avoid alias clash with derived table */
                       select  ind.case_id
                             , ind.sort_id
                             , loi.ind_pref_tera  as ind_label
                       from    case_ind ind
                       inner join lo_indications loi
                           on  loi.ind_pref_tera = ind.ind_pref_tera
                       where   ind.deleted is null
                   ) ind_sub
            on  ind_sub.case_id = ca.case_id
            and ind_sub.sort_id = cp.sort_id

        inner join case_event ce
            on  ca.case_id   = ce.case_id
            and ce.deleted  is null

        /* ── WHERE ──────────────────────────────────────────────────── */
        where
                ca.deleted        is null
            and (ca.state_id <> 1 or cp.sort_id <> 1)  /* Bug 1 – was "(or cp.sort_id-1)" → ORA-00936 */
            and ca.init_rcpt_date  < DATE '2025-07-04'

            /* ── DRUG FILTER: neridronate / neridronato di sodio ─────── */
            and (
                    upper(dbms_lob.substr(cp.product_name, 200, 1))
                        like '%NERIDRONATE SODIUM%'
                 or upper(dbms_lob.substr(cp.generic_name, 200, 1))
                        like '%NERIDRONATE SODIUM%'
                 or upper(dbms_lob.substr(cp.product_name, 200, 1))
                        like '%NERIDRONATO%'
                 or upper(dbms_lob.substr(cp.generic_name, 200, 1))
                        like '%NERIDRONATO%'
                 or upper(dbms_lob.substr(cp.product_name, 200, 1))
                        like '%NERIDRONIC ACID%'
                 or upper(dbms_lob.substr(cp.generic_name, 200, 1))
                        like '%NERIDRONIC ACID%'
            )

            /* ── EVENT FILTER: off-label indication only ─────────────── */
            and lower(ce.pref_tera) = 'off label use'

        order by ca.case_id, ce.sort_id   /* Bug 5 – was ce.seq_num (non-existent col); using sort_id */
    );

    disconnect from oracle;
quit;
