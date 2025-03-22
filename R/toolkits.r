

########加载全局使用包############
usethis::use_package(package = "dplyr", type = "depends")
usethis::use_package(package = "data.table", type = "depends")
usethis::use_package(package = "ggplot2", type = "depends")
usethis::use_package(package = "Rtsne", type = "Imports")
usethis::use_package(package = "stats", type = "depends")
usethis::use_package(package = "ggsci", type = "Imports")
usethis::use_package(package = "edgeR", type = "Imports")
usethis::use_package(package = "limma", type = "Imports")
usethis::use_package(package = "ggpubr", type = "Imports")
#' Title a small function to import dataframe in a exp manner
#'
#' @param dir the dir of data.frame
#' @param col the name of colnames
#' @keywords internal
#' @return a exp
#' @examples
.T_fread_dataframe <- function(dir,col="V1"){
  A <- fread(dir)%>%as.data.frame()
  row.names(A)<- A[,col]
  A[,col] <- NULL
  A <- as.data.frame(A)
  A
}






#' Title .T_ENST2TABLE
#' covert ENST to TABLE
#' @param ENST_res an ENST vector
#'
#' @return table with ENST ENSG genename
#' @keywords internal

.T_ENST2TABLE<- function(ENST_res){
  data(enst_ensg_genename)
  if(!is.vector(ENST_res)){
    stop("input must be vector!")
  }

  ENST_res <- data.table::as.data.table(ENST_res)
  colnames(ENST_res) <- "ENST"
table <- data.table::merge.data.table(ENST_res,enst_ensg_genename,by.x = "ENST",by.y = "ensembl_transcript_id",all.x = T)
table
}







