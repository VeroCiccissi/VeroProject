/* =============================================================================
   ARGUS SAFETY - Estrazione Principio Attivo Unico per Segnalazione
   Scopo   : Restituisce un solo principio attivo (GENERIC_NAME) per ciascuna
             segnalazione, eliminando i duplicati causati dai followup.
   Filtri  : Paracetamolo, DRUG_TYPE=1 (Suspect), GEN-FEB 2026,
             esclude prodotti "Not company"
   DB      : Oracle Argus Safety
   Nota    : La deduplicazione avviene su CASE_ID.
             Il DISTINCT esterno elimina i duplicati da join con CASE_FOLLOWUP.
             Sostituire le date e il nome del principio attivo secondo necessità.
   ============================================================================= */

SELECT DISTINCT
    cm.case_id,
    cm.case_num,
    UPPER(cp.generic_name)   AS principio_attivo,
    cp.product_name          AS nome_prodotto,
    cm.init_rept_date        AS data_ricezione_iniziale
FROM
    case_master    cm,
    case_product   cp,
    case_followup  cf
WHERE
    /* Join principale */
    cm.case_id  = cp.case_id
    /* Join outer con followup (una segnalazione può non avere followup) */
    AND cm.case_id = cf.case_id (+)
    /* Solo farmaco sospetto primario */
    AND cp.drug_type = 1
    /* Principio attivo = Paracetamolo */
    AND UPPER(cp.generic_name) LIKE UPPER('%Paraceta%')
    /* Esclude prodotti non aziendali */
    AND UPPER(cp.product_name) NOT LIKE UPPER('%Not company%')
    /* Filtro date: segnalazione iniziale OPPURE followup nel periodo */
    AND (
            (
                TRUNC(cm.init_rept_date) >= TO_DATE('01-JAN-2026', 'DD-MON-YYYY')
                AND TRUNC(cm.init_rept_date) <= TO_DATE('28-FEB-2026', 'DD-MON-YYYY')
            )
            OR
            (
                TRUNC(cf.receipt_date)   >= TO_DATE('01-JAN-2026', 'DD-MON-YYYY')
                AND TRUNC(cf.receipt_date)   <= TO_DATE('28-FEB-2026', 'DD-MON-YYYY')
            )
       )
ORDER BY
    cm.case_num;
