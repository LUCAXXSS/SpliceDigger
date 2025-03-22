

#' Title SpliceDigger.visualize_AS_transcript
#'
#' @param SD_object
#' @param feature
#' @param feature_type
#' @param anno_to_use
#' @param anno_to_demo
#' @param anno_list
#' @param calculate_method
#' @param color
#'
#' @return
#' @export
#'
#' @examples
SpliceDigger.visualize_AS_transcript <- function(SD_object=A,feature="EGFR",feature_type=c("ENST", "Symbol", "ENSG"),
                                                          anno_to_use="sort",
                                                          anno_to_demo="BCBM",
                                                          anno_list =result,
                                                          calculate_method=c("mean", "median"),
                                                          color = "#00A087"){
  ##pheno是两列的矩阵，一列是样本名，一列是cluster名
  ##results是列名对应的ENST

  feature_type=match.arg(feature_type,c("ENST", "Symbol", "ENSG"))
  calculate_method=match.arg(calculate_method,c("mean", "median"))
  pheno <- SD_object$anno[anno_to_use]
  pheno_use <- pheno[pheno[[anno_to_use]]%in%anno_to_demo,,drop=F]


  exp <- SD_object$iso_exp_raw

  exp_use <- exp[,colnames(exp)%in%row.names(pheno_use)]




  data("enst_ensg_genename")


  #return standard ENSG ID from three sources
  if(feature_type=="ENST"){
    input_ENSG <- enst_ensg_genename$ensembl_gene_id[enst_ensg_genename$ensembl_transcript_id==feature]%>%unique()
    if(!length(input_ENSG)>0){
      stop("ENST ID not found in database!")
    }

  }

  if(feature_type=="Symbol"){
    input_ENSG <- enst_ensg_genename$ensembl_gene_id[enst_ensg_genename$external_gene_name==feature]%>%unique()
    if(!length(input_ENSG)>0){
      stop("Symbol not found in database!")
    }
  }

  if(feature_type=="ENSG"){
    input_ENSG <- enst_ensg_genename$ensembl_gene_id[enst_ensg_genename$external_gene_name==feature]%>%unique()
    if(!length(input_ENSG)>0){
      stop("ENSG ID not found in database!")
    }
  }


  ENST.use <- enst_ensg_genename$ensembl_transcript_id[enst_ensg_genename$ensembl_gene_id==input_ENSG]
  Symbol.use <- enst_ensg_genename$external_gene_name[enst_ensg_genename$ensembl_gene_id==input_ENSG]%>%unique()
  exp_use <- exp_use[row.names(exp_use)%in%ENST.use,]%>%t()%>%as.data.frame()

  #如果一个基因完全没有表达（ENSG上也没有表达的话），其比例无法计算，所以设定为负值，在统计中改为0#
  exp_use[exp_use<0] <- 0

  if (calculate_method=="mean") {
  exp_use_plot <- data.frame(colnames(exp_use),colMeans(exp_use))
  }
  if (calculate_method=="median") {
    exp_use_plot <- data.frame(colnames(exp_use),sapply(exp_use, median))
  }
  colnames(exp_use_plot) <- c("ID","Value")

  ###plot module#############

  p <- ggpubr::ggbarplot(exp_use_plot,
                 x = "ID",
                 y = "Value",
                 fill = color,
                 color = "black",
                 palette = NULL,title = paste0(" ALL ENST in ",input_ENSG," GeneSymbol:",Symbol.use)
  )+ylim(0,1)+coord_flip() # 禁用默认调色板



  p
}




