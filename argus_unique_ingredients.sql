/* =============================================================================
   ARGUS SAFETY - Estrazione Principio Attivo Unico per Segnalazione
   Scopo : Restituisce un solo principio attivo (ingrediente primario) per
           ciascuna segnalazione (case), eliminando i duplicati.
   DB    : Oracle Argus Safety
   Nota  : Nel caso di più principi attivi per la stessa segnalazione viene
           mantenuto quello con il rank più basso (SEQ_NUM minimo del prodotto
           sospettato primario). Modificare la clausola ORDER BY per adottare
           un criterio diverso.
   ============================================================================= */

SELECT
    case_num,
    case_id,
    report_type,
    report_date,
    product_name,
    ingredient_name,
    ingredient_code
FROM (
    SELECT
        cm.case_num,
        cm.case_id,
        cm.report_type,
        cm.date_receiv                              AS report_date,
        lp.product_name,
        UPPER(TRIM(li.ingredient_name))             AS ingredient_name,
        li.ingredient_id                            AS ingredient_code,
        ROW_NUMBER() OVER (
            PARTITION BY cm.case_id
            ORDER BY
                /* Prodotto sospettato primario (drug_type = 1) prima degli altri */
                CASE WHEN cp.drug_type = 1 THEN 0 ELSE 1 END,
                /* A parità di tipo, usa il numero di sequenza del prodotto */
                cp.seq_num,
                /* A parità assoluta, ordina alfabeticamente sul nome ingrediente */
                UPPER(TRIM(li.ingredient_name))
        )                                           AS rn
    FROM
        case_master          cm
        /* Prodotti associati alla segnalazione */
        JOIN case_product     cp  ON cp.case_id       = cm.case_id
        /* Anagrafica prodotto */
        JOIN lm_product       lp  ON lp.product_id    = cp.product_id
        /* Tabella di raccordo prodotto-ingrediente */
        JOIN lm_prod_ingred   pi  ON pi.product_id    = lp.product_id
        /* Anagrafica ingrediente / principio attivo */
        JOIN lm_ingredients   li  ON li.ingredient_id = pi.ingredient_id
    WHERE
        /* Solo casi validi (non cancellati) */
        cm.case_deleted  = 0
        /* Solo prodotti con almeno un ruolo di farmaco (esclude dispositivi, ecc.) */
        AND cp.drug_type IN (1, 2, 3)   -- 1=Suspect, 2=Concomitant, 3=Interacting
        /* Solo ingredienti attivi in anagrafica */
        AND li.active_flag = 1
)
WHERE rn = 1   -- conserva solo il primo principio attivo per segnalazione
ORDER BY
    case_num;
