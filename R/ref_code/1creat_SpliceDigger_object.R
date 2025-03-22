

#' Title
#'
#' @param ENST_table
#' @param exp_table
#' @param calculate_pct
#'
#' @return
#' @export
#'
#' @examples
SpliceDigger.create_SpliceDigger_object <- function(ENST_table=iso_tpm_You_etal_2025,exp_table=exp_You_etal_2025,calculate_pct=T){
  if (!is.data.frame(ENST_table)) {
    stop(paste0(substitute(ENST_table)," must be a data.frame"))
  }

  if (!is.data.frame(exp_table)) {
    stop(paste0(substitute(exp_table)," must be a data.frame"))
  }

  if (any(duplicated(colnames(ENST_table)))) {
    stop(paste0(substitute(ENST_table), " colnames not unique!"))
  }

  if (any(duplicated(colnames(exp_table)))) {
    stop(paste0(substitute(exp_table), " colnames not unique!"))
  }

  ####create SpliceDiggerObject

  ENST_table <- as.data.frame(ENST_table)
  exp_table <- as.data.frame(exp_table)

  SD.list <- list(iso_exp_raw=ENST_table,exp_raw=exp_table)
  class(SD.list) <- "SpliceDiggerObject"

  #######calculate pct_table########


  if (calculate_pct) {


    iso_exp_pct_raw <- SD.list[["iso_exp_raw"]]
    iso_exp_pct_raw$ENST <- row.names(iso_exp_pct_raw)


    ENST_list <- SpliceDigger:::.T_ENST2TABLE(row.names(iso_exp_pct_raw))

    iso_exp_pct_calculate <- merge(ENST_list,iso_exp_pct_raw,by = "ENST",all.y = T)


    iso_exp_pct <- iso_exp_pct_calculate %>%
      group_by(ensembl_gene_id) %>%
      mutate(across(where(is.numeric), ~ . / sum(.), .names = "{.col}"))%>%as.data.frame()



    ######数据格式整理清洗
    row.names(iso_exp_pct) <- iso_exp_pct$ENST

    iso_exp_pct <- iso_exp_pct[,!colnames(iso_exp_pct)%in%colnames(ENST_list)]
    iso_exp_pct <- as.matrix(iso_exp_pct)


    iso_exp_pct[is.nan(iso_exp_pct)] <- NA
    iso_exp_pct <- as.data.frame(iso_exp_pct)

    SD.list[["iso_exp_pct"]] <- iso_exp_pct
    SD.list
  }







}




