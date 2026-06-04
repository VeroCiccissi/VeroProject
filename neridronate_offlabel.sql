/* ============================================================================
   neridronate_offlabel.sql
   ------------------------
   Estrazione casi AE (NERIDRONATE – uso off-label) dal database di
   farmacovigilanza Oracle PROD1 (PISAS01.abiogen.it) via SAS passthrough.

   Tabella output : NERIDRONATE_OFFLABEL
   Data lock      : 2025-07-04
   Step           : STEP 3 N – VERSIONE FINALE (FIX CLOB + ORA-00936)
   ============================================================================ */

proc sql;
    connect to oracle (user="&user" password="&pw" path="PROD1!");

    create table NERIDRONATE_OFFLABEL as
    select *
    from connection to oracle
    (
        select distinct
              ca.case_id
            , ca.case_seq_num
            , ca.init_rcpt_date                         as init_receipt_date
            , dbms_lob.substr(cp.product_name, 200, 1)  as cp_product_name
            , cpi.pat_age
            , cpi.gender_id
            , lgo.group_name                            as age_group
            , cp.sort_id
            , cp.pref_tera
            , dbms_lob.substr(cp.generic_name, 200, 1)  as generic_name
            , case
                  when lower(ce.pref_tera) = 'off label use'
                  then ce.pref_tera
              end                                       as PT
            , ce.body_sys                               as soc

        from case_master ca

        inner join case_pat_info cpi
            on  ca.case_id   = cpi.case_id
            and cpi.deleted  is null

        left  join lo_age_groups lgo
            on  lgo.age_group_id = cpi.age_group_id

        left  join lo_gender lg
            on  lg.gender_id     = cpi.gender_id

        inner join case_product cp
            on  ca.case_id      = cp.case_id
            and cp.deleted      is null
            and (cp.first_tera  is null or cp.first_tera <> 2)

        left  join (
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
            on  ca.case_id  = ce.case_id
            and ce.deleted  is null

        where
                ca.deleted         is null
            and ca.state_id        <> 1
            and cp.sort_id         <> 1
            and ca.init_rcpt_date   < DATE '2025-07-04'

            /* DRUG FILTER – neridronato di sodio */
            and (
                    upper(dbms_lob.substr(cp.product_name, 200, 1)) like '%NERIDRONATE SODIUM%'
                 or upper(dbms_lob.substr(cp.generic_name, 200, 1)) like '%NERIDRONATE SODIUM%'
                 or upper(dbms_lob.substr(cp.product_name, 200, 1)) like '%NERIDRONATO%'
                 or upper(dbms_lob.substr(cp.generic_name, 200, 1)) like '%NERIDRONATO%'
                 or upper(dbms_lob.substr(cp.product_name, 200, 1)) like '%NERIDRONIC ACID%'
                 or upper(dbms_lob.substr(cp.generic_name, 200, 1)) like '%NERIDRONIC ACID%'
            )

            /* EVENT FILTER – solo indicazione off-label */
            and lower(ce.pref_tera) = 'off label use'

        order by ca.case_id, ce.sort_id
    );

    disconnect from oracle;
quit;
