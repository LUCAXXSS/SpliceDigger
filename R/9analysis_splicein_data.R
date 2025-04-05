

#' Title
#'
#' @param SD_object
#' @param return_pie_plot
#' @param anno_to_use
#' @param cluster_to_check
#'
#' @return
#' @export
#'
#' @examples
SpliceDigger.plot_AS_pieplot <- function(SD_object=A,return_pie_plot=T,
                                         anno_to_use="iso_exp_pct_filtered_kmeans_cluster",cluster_to_check="cluster_2"){

sample_Select <- row.names(SD_object$anno)[SD_object$anno[[anno_to_use]]==cluster_to_check]

 psi_exp <- SD_object$psi[["psi_matrix"]]
 psi_anno<- SD_object$psi[["psi_anno"]]

 psi_exp_select_RAW <- psi_exp[,colnames(psi_exp)%in%sample_Select]



 pct_AS_forplot <- data.frame()
 for (i in unique(psi_anno$type)) {
 psi_anno_select <- psi_anno[psi_anno$type==i,]
 psi_exp_select <- psi_exp_select_RAW[psi_anno$type==i,]%>%as.matrix()

 psi_exp_select[!is.numeric(psi_exp_select)] <- NA
 psi_exp_select[is.nan(psi_exp_select)] <- NA

 psi_value<- median(psi_exp_select,na.rm = T)  ####使用中位数

 B<- data.frame(type=i,value=psi_value)
 pct_AS_forplot <- rbind(pct_AS_forplot,B)

 }
 pct_AS_forplot <-pct_AS_forplot[!is.nan(pct_AS_forplot$value),]
 pct_AS_forplot <- pct_AS_forplot %>%
   mutate(
     # 计算累计百分比用于标签位置
     cumsum = cumsum(value),
     pos = cumsum - value/2,
     # 创建百分比标签文本
     label = paste0(round(value/sum(value)*100, 1), "%")
   )

  p <- ggplot2::ggplot(pct_AS_forplot, aes(x="type", y=value, fill=type)) +
   geom_bar(stat="identity", width=1, color="white") +
   coord_polar("y", start=0) +
   theme_void()+ # remove background, grid, numeric labels
    ggrepel::geom_text_repel(aes(y = pos, label = label),
                    color = "white",
                    point.padding = NA)+ggsci::scale_fill_lancet()



  ifelse(return_pie_plot,return(p),return(pct_AS_forplot))

}




#' Title
#'
#' @param SD_object
#' @param return_pie_plot
#' @param gene_to_use
#' @param anno_to_use
#' @param cluster_to_check
#'
#' @return
#' @export
#'
#' @examples
SpliceDigger.plot_gene_AS_pieplot <- function(SD_object=A,return_pie_plot=T,gene_to_use=c("EGFR"),
                                         anno_to_use="iso_exp_pct_filtered_kmeans_cluster",cluster_to_check="cluster_2"){

  sample_Select <- row.names(SD_object$anno)[SD_object$anno[[anno_to_use]]==cluster_to_check]

  psi_exp <- SD_object$psi[["psi_matrix"]]
  psi_anno<- SD_object$psi[["psi_anno"]]


  psi_row_to_keep <-  which(psi_anno$external_gene_name %in% gene_to_use)



  psi_exp <- psi_exp[psi_row_to_keep,]
  psi_anno <- psi_anno[psi_row_to_keep,]

  psi_exp_select_RAW <- psi_exp[,colnames(psi_exp)%in%sample_Select]



  pct_AS_forplot <- data.frame()
  for (i in unique(psi_anno$type)) {
    psi_anno_select <- psi_anno[psi_anno$type==i,]
    psi_exp_select <- psi_exp_select_RAW[psi_anno$type==i,]%>%as.matrix()

    psi_exp_select[!is.numeric(psi_exp_select)] <- NA
    psi_exp_select[is.nan(psi_exp_select)] <- NA

    psi_value<- mean(psi_exp_select,na.rm = T)  ####使用平均数

    B<- data.frame(type=i,value=psi_value)
    pct_AS_forplot <- rbind(pct_AS_forplot,B)

  }

  pct_AS_forplot <-pct_AS_forplot[!is.nan(pct_AS_forplot$value),]

  pct_AS_forplot <- pct_AS_forplot %>%
    mutate(
      # 计算累计百分比用于标签位置
      cumsum = cumsum(value),
      pos = cumsum - value/2,
      # 创建百分比标签文本
      label = paste0(round(value/sum(value)*100, 1), "%")
    )



    p <- ggplot2::ggplot(pct_AS_forplot, aes(x="type", y=value, fill=type)) +
      geom_bar(stat="identity", width=1, color="white") +
      coord_polar("y", start=0) +
      theme_void()+ # remove background, grid, numeric labels
      ggrepel::geom_text_repel(aes(y = pos, label = label),
                               color = "white",
                               point.padding = NA)+ggsci::scale_fill_lancet()

   ifelse(return_pie_plot,return(p),return(pct_AS_forplot))





}
