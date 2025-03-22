

#' Title SpliceDigger.filter_SDobject
#'  scale and filter HVGs
#' @param SD_object
#' @param AS_threshold
#' @param AS_var_topn
#' @param exp_threshold
#' @param exp_var_topn
#'
#' @return
#' @export
#'
#' @examples
SpliceDigger.filter_SDobject<- function(SD_object=A,AS_threshold=0.1,AS_var_topn=5000,
                                                  exp_threshold=20,exp_var_topn=5000){

  iso_exp <- SD_object$iso_exp_pct
  iso_exp[is.na(iso_exp)] <- 0

  ##########Filter HVGs_AS#########
  iso_exp <- iso_exp[rowMeans(iso_exp)>AS_threshold,]
  # 第三步：生成一个方差表，并确保行名保留


  exp_var <- apply(iso_exp, 1, var)  # 直接计算方差并保留为一个向量
  exp_var <- data.frame(var = exp_var, row.names = names(exp_var))  # 将其转换为数据框，行名重建
  exp_var$keep <- "keep"
  exp_var <- exp_var[order(exp_var$var, decreasing = TRUE),]
  exp_var$keep <-NULL

  if (length(row.names(exp_var))<exp_var_topn) {
    stop("AS_threshold or AS_var_topn too high,please adjust one of them!")
  }


  gene_keep <- row.names(exp_var)[1:AS_var_topn]

  iso_exp <- iso_exp[gene_keep,]

  ##########Filter HVGs_exp#########
  exp <- SD_object$exp_raw
  exp[is.na(exp)] <- 0

  exp <- exp[rowMeans(exp)>exp_threshold,]
  # 第三步：生成一个方差表，并确保行名保留
  exp_var <- apply(exp, 1, var)  # 直接计算方差并保留为一个向量
  exp_var <- data.frame(var = exp_var, row.names = names(exp_var))  # 将其转换为数据框，行名重建
  exp_var$keep <- "keep"
  exp_var <- exp_var[order(exp_var$var, decreasing = TRUE),]
  exp_var$keep <-NULL

  if (length(row.names(exp))<AS_var_topn) {
    stop("exp_threshold or exp_var_topn too high,please agjust one of them!")
  }


  gene_keep <- row.names(exp_var)[1:exp_var_topn]

  exp <- exp[gene_keep,]

  SD_object$iso_exp_pct_filtered <- iso_exp
  SD_object$exp_filtered <- exp
  SD_object
}





