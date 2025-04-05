
#' Title SpliceDigger.add_SUPPA2_psi
#'
#' @param SD_object
#' @param psi_matrix
#'
#' @return
#' @export
#'
#' @examples
SpliceDigger.add_SUPPA2_psi <- function(SD_object,psi_matrix){

  if (!class(SD_object)=="SpliceDiggerObject") {
    stop(paste0(substitute(SD_object)," must be a SpliceDiggerObject"))
  }

  nameorder <- colnames(SD_object[["iso_exp_raw"]])
  if (!identical(colnames(psi_matrix),nameorder)) {
    stop(paste0(substitute(SD_object)," and ",substitute(psi_matrix)," must be have identical colnames!",
                "\n current ",substitute(SD_object)," colnames are ",toString(nameorder)))

  }

  #save data to SDobject

  SD_object$psi[["psi_matrix"]] <- psi_matrix

  #create psi anno from colnames of psi
  psi_anno_raw <- data.frame(row.names(psi_matrix))



  # extract data from row.names
  split_result <- strsplit(psi_anno_raw[,1], ";")

  # 使用 lapply 和 sapply 批量提取 feature1 和 feature2
  feature1 <- sapply(split_result, function(x) x[1])                             # 提取第一个字段
  feature2 <- sapply(split_result, function(x) strsplit(x[2], ":")[[1]][1])     # 提取第二个字段

  # 构造数据框
  AS_events_anno <- data.frame(ID=psi_anno_raw[,1],ENSG = feature1, type = feature2, stringsAsFactors = FALSE)
  data("enst_ensg_genename")

  AS_events_anno <- merge(AS_events_anno,enst_ensg_genename,by.x = "ENSG",by.y="ensembl_gene_id",all.x = T,no.dups = T)
  AS_events_anno <- AS_events_anno[!duplicated(AS_events_anno$ID),]

  AS_events_anno$ID <- NULL
  AS_events_anno$ensembl_transcript_id <- NULL

  SD_object$psi[["psi_anno"]] <- AS_events_anno
  SD_object
}
