/* ============================================================================
   neridronate_offlabel.sql
   ------------------------
   AE_CASES extraction for NERIDRONATE (off-label use).
   SAS/ACCESS Oracle passthrough  |  PISAS01.abiogen.it  |  DB: PROD1
   Output table : NERIDRONATE_OFFLABEL
   Data lock    : 2025-07-04

   STEP 3 N – FINAL VERSION (FIX CLOB ERROR + ORA-00936)

   Corrections vs original
   -----------------------
   1. ORA-00936 – "(or cp.sort_id-1)" was not valid SQL; corrected to
      "(ca.state_id <> 1 or cp.sort_id <> 1)".
   2. dbms_lob.substr – explicit 3-argument form (lob, amount, offset => 1)
      used throughout; 2-arg form is ambiguous in Oracle SQL passthrough.
   3. ca.cp_product_name – column does not exist on CASE_MASTER; replaced
      with dbms_lob.substr(cp.product_name, 200, 1) from CASE_PRODUCT.
   4. CASE ... THEN cp.pref_tera – the PT column should carry the event
      preferred term, not the product term; corrected to ce.pref_tera.
   5. ind_pref_tera subquery – inner column renamed to ind_label to avoid
      an alias clash with the derived-table identifier.
   6. ORDER BY ce.seq_num – column does not exist on CASE_EVENT;
      corrected to ce.sort_id.
   ============================================================================ */

proc sql;
    connect to oracle (user="&user" password="&pw" path="PROD1!");

    /* ── CREATE OUTPUT TABLE ─────────────────────────────────────────────── */
    create table NERIDRONATE_OFFLABEL as
    select *
    from connection to oracle
    (
        select distinct
              ca.case_id
            , ca.case_seq_num
            , ca.init_rcpt_date                          as init_receipt_date
            /* cp.product_name is a CLOB – read via dbms_lob.substr         */
            , dbms_lob.substr(cp.product_name,  200, 1)  as cp_product_name
            , cpi.pat_age
            , cpi.gender_id
            , lgo.group_name                             as age_group
            , cp.sort_id
            , cp.pref_tera
            /* generic_name is a CLOB – read via dbms_lob.substr            */
            , dbms_lob.substr(cp.generic_name,  200, 1)  as generic_name
            /* PT: event preferred term when the event is an off-label use  */
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
            on  ca.case_id       = cp.case_id
            and cp.deleted       is null
            and (cp.first_tera   is null or cp.first_tera <> 2)

        left  join (
                       /* preferred indication term matched to product-sort  */
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
            and ce.deleted   is null

        /* ── FILTERS ────────────────────────────────────────────────────── */
        where
                ca.deleted        is null
            and (ca.state_id <> 1 or cp.sort_id <> 1)   /* fix 1: was "(or cp.sort_id-1)" */
            and ca.init_rcpt_date  < DATE '2025-07-04'

            /* DRUG FILTER – neridronate / neridronato di sodio             */
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

            /* EVENT FILTER – off-label use indication only                 */
            and lower(ce.pref_tera) = 'off label use'

        order by ca.case_id, ce.sort_id   /* fix 6: was ce.seq_num */
    );

    disconnect from oracle;
quit;
