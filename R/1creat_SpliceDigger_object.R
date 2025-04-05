

#' Title SpliceDigger.create_SD_object
#' return a list object with ENST_raw,exp_raw and ENSG_pct
#' @param ENST_table A data.frame where row.names be ENST names and colnames be samplenames, this
#'         file could be usually generated from SUPPA2 pipelines
#' @param exp_table A data.frame where row.names be gene names and colnames be samplenames, this
#'         file could be usually generated from expression matrix pipelines
#' @param calculate_pct if calculate all ENST percentage under one ENSG number
#'
#' @return A SpliceDigger object
#' @export
#'
#' @examples
#'

SpliceDigger.create_SD_object <- function(ENST_table=iso_tpm_You_etal_2025,exp_table=exp_You_etal_2025,calculate_pct=T){
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



#' Title SpliceDigger.add_annotation
#' Add annotation to SpliceDigger object
#' @param SD_object
#' @param anno a data.frame containing all annotaion info, row.names should be set to sample names
#'
#' @return
#' @export
#'
#' @examples
SpliceDigger.add_annotation <- function(SD_object,anno){

  if (!is.data.frame(anno)) {
    stop(paste0(substitute(anno)," must be a data.frame"))
  }

  nameorder <- colnames(SD_object[["iso_exp_raw"]])
  if (!identical(row.names(anno),nameorder)) {
    stop(paste0(substitute(SD_object)," and ",substitute(anno)," must be have identical colnames!",
                "\n current ",substitute(SD_object)," colnames are ",toString(nameorder)))

  }

  if(length(SD_object[["anno"]])>0){
    colname<-colnames(SD_object[["anno"]])
    colname_merge <- c(colname,colnames(anno))
    SD_object[["anno"]] <-cbind(SD_object[["anno"]],anno)
    colnames(SD_object[["anno"]]) <- colname_merge
    SD_object[["anno"]] <- SD_object[["anno"]][,!duplicated(colnames(SD_object[["anno"]])),drop = FALSE]
  }else{
    SD_object[["anno"]] <-anno

  }
  SD_object

}


