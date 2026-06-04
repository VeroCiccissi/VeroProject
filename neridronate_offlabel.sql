/* ============================================================================
   neridronate_offlabel.sql
   ------------------------
   Adverse-event case extraction for NERIDRONATE (off-label use).
   Run via SAS/ACCESS passthrough against the Oracle PROD1 safety DB
   (PISAS01.abiogen.it).

   Step : STEP 3 N – FINAL VERSION (FIX CLOB ERROR)
   Output table : NERIDRONATE_OFFLABEL
   Data lock     : 2025-07-04
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
            , ca.cp_product_name
            , ca.case_seq_num
            , ca.init_rcpt_date                        as init_receipt_date
            /* cp.product_name  */
            , cpi.pat_age
            , cpi.gender_id
            , lgo.group_name                           as age_group
            , cp.sort_id
            , cp.pref_tera

            /* CLOB columns read via dbms_lob.substr to avoid CLOB errors */
            , dbms_lob.substr(cp.generic_name, 200)    as generic_name

            /* A EVENT (OFF LABEL) */
            , case
                  when lower(ce.pref_tera) = 'off label use'
                  then cp.pref_tera
              end                                      as PT

            , ce.body_sys                              as soc

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
                       /* Preferred indication term per product-sort */
                       select  ind.case_id
                             , ind.sort_id
                             , loi.ind_pref_tera
                       from    case_ind ind
                       inner join lo_indications loi
                           on  loi.ind_pref_tera = ind.ind_pref_tera
                       where   ind.deleted is null
                   ) ind_pref_tera
            on  ind_pref_tera.case_id = ca.case_id
            and ind_pref_tera.sort_id = cp.sort_id

        inner join case_event ce
            on  ca.case_id   = ce.case_id
            and ce.deleted  is null

        /* ── WHERE ──────────────────────────────────────────────────── */
        where
                ca.deleted       is null
            and ca.state_id      <> 1
            and cp.sort_id       <> 1          /* fixed: was "(or cp.sort_id-1)" → ORA-00936 */
            and ca.init_rcpt_date < DATE '2025-07-04'

            /* ── DRUG FILTER: neridronate / neridronato di sodio ─────── */
            and (
                    upper(dbms_lob.substr(cp.product_name, 200))
                        like '%NERIDRONATE SODIUM%'
                 or upper(dbms_lob.substr(cp.generic_name, 200))
                        like '%NERIDRONATE SODIUM%'
                 or upper(dbms_lob.substr(cp.product_name, 200))
                        like '%NERIDRONATO%'
                 or upper(dbms_lob.substr(cp.generic_name, 200))
                        like '%NERIDRONATO%'
                 or upper(dbms_lob.substr(cp.product_name, 200))
                        like '%NERIDRONIC ACID%'
                 or upper(dbms_lob.substr(cp.generic_name, 200))
                        like '%NERIDRONIC ACID%'
            )

            /* ── EVENT FILTER: off-label indication only ─────────────── */
            and lower(ce.pref_tera) = 'off label use'

        order by ca.case_id, ce.seq_num
    );

    disconnect from oracle;
quit;
