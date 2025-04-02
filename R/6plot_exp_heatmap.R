
#' Title SpliceDigger.visualize_exp_heatmap
#'
#' @param SD_object
#' @param use_data
#' @param feature_to_plot a vector of gene names, typically c("SRSF1","SRSF2")
#' @param feature_type
#' @param anno_to_use anno to use in SD_object$anno
#' @param FDR_threshold default1
#' @param color a vector, default c("navy", "white", "firebrick3")
#'
#' @return
#' @export
#'
#' @examples
SpliceDigger.visualize_exp_heatmap <- function(SD_object,use_data="exp_raw",feature_to_plot,
                                                 anno_to_use,FDR_threshold=1,
                                                 color = c("navy", "white", "firebrick3")){

  exp <- SD_object$exp_raw
  genenames <- row.names(exp)%>%as.data.frame()
  colnames(genenames) <- "symbol"


  Var_deg_res <- list()
  for (i in unique(SD_object$anno[[anno_to_use]])) {
  df <- SpliceDigger.findMarkers(SD_object =SD_object, anno_to_use = anno_to_use,use_data = use_data,cluster_to_check = i,only.pos = F,FDR.threshold = FDR_threshold,logfc.threshold = 0)
  res <- df["logFC"]
  res$symbol <- row.names(res)
  res_all <- genenames
  res_all <- merge(res_all,res,by = "symbol",all.x = T)

  res_small <-data.frame( res_all["logFC"],row.names = res_all$symbol)
  colnames(res_small) <- i
  Var_deg_res[[i]] <- res_small
  }

  exp_logfc <- as.data.frame(Var_deg_res)


  exp_plot <- exp_logfc[row.names(exp_logfc)%in%feature_to_plot,]
  exp_plot[is.na(exp_plot)] <- 0

  p <- pheatmap::pheatmap(exp_plot,scale = "none",cluster_rows = F,cluster_cols = F,treeheight_col =0,fontsize_row = 5,
           color = colorRampPalette(color)(50))

  p
}
