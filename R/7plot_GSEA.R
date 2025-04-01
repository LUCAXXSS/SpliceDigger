
#' Title
#'
#' @param SD_object
#' @param use_data
#' @param anno_to_use
#' @param cluster_to_check
#' @param topn
#' @param msigdbr_species
#' @param msigdbr_category
#' @param msigdbr_subcategory
#' @param gsea_pvalue
#' @param plot
#' @param color
#'
#' @return
#' @export
#'
#' @examples
SpliceDigger.plot_GSEA <- function(SD_object=A,use_data="exp_raw",
                                     anno_to_use="iso_exp_pct_filtered_kmeans_cluster",cluster_to_check="cluster_2",
                                   topn=20,
                                   msigdbr_species = "Homo sapiens",
                                   msigdbr_category="H",msigdbr_subcategory=NULL,gsea_pvalue=0.05,plot=T,color="#DC0000"){


  if (!length(msigdbr_subcategory>0)) {

    gene_sets <- msigdbr::msigdbr(species = msigdbr_species,category = msigdbr_category)
  }else{

  gene_sets <- msigdbr::msigdbr(species = msigdbr_species,category = msigdbr_category, subcategory = msigdbr_subcategory)

}

  gene_sets <- select(gene_sets,c("gs_name", "gene_symbol"))


  df <- SpliceDigger.findMarkers(SD_object,use_data =use_data ,anno_to_use = anno_to_use,cluster_to_check = cluster_to_check,
                                 only.pos = F)
  print(length(row.names(df)))

  gene_to_use <- df$logFC
  names(gene_to_use) <- row.names(df)
  gene_to_use
  gsea_res <- clusterProfiler::GSEA(gene_to_use,
                   TERM2GENE = gene_sets,
                   minGSSize = 1,
                   maxGSSize = 500,
                   pvalueCutoff = gsea_pvalue,
                   pAdjustMethod = "BH",
                   seed = 123
  )

 gsea_res_forplot <- gsea_res@result
 gsea_res_forplot <- gsea_res_forplot[order(gsea_res_forplot$NES,decreasing = T),]

 gsea_res_forplot_final <- gsea_res_forplot[topn:1,]

 gsea_res_forplot_final

 if(plot){
   gsea_res_forplot_final <- na.omit(gsea_res_forplot_final)

   gsea_res_forplot_final$p.value <- gsea_res_forplot_final$p.adjust


   p <- ggpubr::ggbarplot(gsea_res_forplot_final,
                          x = "ID",
                          y = "NES",
                          fill = "p.value",
                          color = "black",
                          palette = NULL
   )+  coord_flip()+scale_fill_gradient(low = color,high = "grey90")

   p

 }



}





SpliceDigger.plot_GSEA_ridgeplot <- function(SD_object=A,use_data="exp_raw",
                                   anno_to_use="iso_exp_pct_filtered_kmeans_cluster",cluster_to_check="cluster_2",
                                   topn=20,custom=F,
                                   msigdbr_species = "Homo sapiens",
                                   msigdbr_category="H",msigdbr_subcategory=NULL,
                                   geneset_to_plot="HALLMARK_ANGIOGENESIS"){


  if (!length(msigdbr_subcategory>0)) {

    gene_sets <- msigdbr::msigdbr(species = msigdbr_species,category = msigdbr_category)
  }else{

    gene_sets <- msigdbr::msigdbr(species = msigdbr_species,category = msigdbr_category, subcategory = msigdbr_subcategory)

  }

  gene_sets <- select(gene_sets,c("gs_name", "gene_symbol"))

  gene_sets <-gene_sets[gene_sets$gs_name==geneset_to_plot,]


  df <- SpliceDigger.findMarkers(SD_object,use_data =use_data ,anno_to_use = anno_to_use,cluster_to_check = cluster_to_check,
                                 only.pos = F)
  print(length(row.names(df)))

  gene_to_use <- df$logFC
  names(gene_to_use) <- row.names(df)
  gene_to_use

  gsea_res <- clusterProfiler::GSEA(gene_to_use,
                                    TERM2GENE = gene_sets,
                                    minGSSize = 1,
                                    maxGSSize = 500,
                                    pvalueCutoff = 1,
                                    pAdjustMethod = "BH",
                                    seed = 123
  )

  gseap1 <- enrichplot::gseaplot2(gsea_res,
                                  geneset_to_plot,#富集的ID编号
                      color = "red", #GSEA线条颜色
                      base_size = 15,#基础字体大小
                      rel_heights = c(1.5, 0.5, 1),#副图的相对高度
                      subplots = 1:3,   #要显示哪些副图 如subplots=c(1,3) #只要第一和第三个图
                      ES_geom = "line", #enrichment score用线还是用点"dot"
                      pvalue_table = T) #显示pvalue等信息

  gseap1



}


