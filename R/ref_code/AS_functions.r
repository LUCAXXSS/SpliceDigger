





#AS_Events_Function 20250306
#Author: Muyuan_You

#Main functions
library(tidyverse)
library(Rtsne) # 载入Rtsne包
library(ggsci)
library(biomaRt)
library(data.table)
library(edgeR)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(biomaRt)

M_filter_highexp_highvar_exp <- function(exp_dataframe,rowmean_threshold=0.1,vartopn=5000){
  exp <- exp_dataframe
  exp <- na.omit(exp)

  ##########筛选高变基因#########
  exp <- exp[rowMeans(exp)>rowmean_threshold,]
  # 第三步：生成一个方差表，并确保行名保留
  exp_var <- apply(exp, 1, var)  # 直接计算方差并保留为一个向量
  exp_var <- data.frame(var = exp_var, row.names = names(exp_var))  # 将其转换为数据框，行名重建
  exp_var$keep <- "keep"
  exp_var <- exp_var[order(exp_var$var, decreasing = TRUE),]
  exp_var$keep <-NULL
  gene_keep <- row.names(exp_var)[1:vartopn]

  exp <- exp[gene_keep,]
  exp
}


M_find_varible_genes <- function(gene_exp,pheno,cluster_to_check,top_n=300){

  ################差异基因分析##############
  exp_deg <- gene_exp
  group = pheno%>%as.vector()


  group[!group==cluster_to_check] <- "ctrl"
  group[group==cluster_to_check] <- "treat"

  group <- factor(group,levels = c("ctrl","treat"))

  dge <- DGEList(counts = exp_deg, group = group)
  # 重新计算库大小
  dge <- calcNormFactors(dge)

  # 估计共同离散度和趋势离散度
  dge <- estimateDisp(dge)

  design.matrix<-model.matrix(~0+group)
  colnames(design.matrix)<-levels(group)

  fit <- glmQLFit(y = dge, design = design.matrix)

  con <- makeContrasts(treat - ctrl, levels=design.matrix)

  #进行差异分析
  qlf <- glmQLFTest(fit, contrast=con)

  dge_results <- as.data.frame(topTags(qlf,n=nrow(dge$counts),sort.by = "logFC"))
  dge_results <- dge_results[order(dge_results$logFC,decreasing = T),]


  dge_results <- dge_results[dge_results$FDR<0.05,]



  GENE_res<- row.names(dge_results)[1:top_n]

  GENE_res


}


M_give_varible_genes_table <- function(gene_exp,pheno,cluster_to_check,top_n=300){

  ################差异基因分析
  exp_deg <- gene_exp
  group = pheno%>%as.vector()


  group[!group==cluster_to_check] <- "ctrl"
  group[group==cluster_to_check] <- "treat"

  group <- factor(group,levels = c("ctrl","treat"))

  dge <- DGEList(counts = exp_deg, group = group)
  # 重新计算库大小
  dge <- calcNormFactors(dge)

  # 估计共同离散度和趋势离散度
  dge <- estimateDisp(dge)

  design.matrix<-model.matrix(~0+group)
  colnames(design.matrix)<-levels(group)

  fit <- glmQLFit(y = dge, design = design.matrix)

  con <- makeContrasts(treat - ctrl, levels=design.matrix)

  #进行差异分析
  qlf <- glmQLFTest(fit, contrast=con)

  dge_results <- as.data.frame(topTags(qlf,n=nrow(dge$counts),sort.by = "logFC"))
  dge_results <- dge_results[order(dge_results$logFC,decreasing = T),]

  dge_results
}


M_find_varible_AS_events <- function(iso_exp,pheno,cluster_to_check,top_n=300){

  ################差异基因分析##############
  exp_deg <- iso_exp
  group = pheno%>%as.vector()


  group[!group==cluster_to_check] <- "ctrl"
  group[group==cluster_to_check] <- "treat"

  group <- factor(group,levels = c("ctrl","treat"))

  dge <- DGEList(counts = exp_deg, group = group)
  # 重新计算库大小
  dge <- calcNormFactors(dge)

  # 估计共同离散度和趋势离散度
  dge <- estimateDisp(dge)

  design.matrix<-model.matrix(~0+group)
  colnames(design.matrix)<-levels(group)

  fit <- glmQLFit(y = dge, design = design.matrix)

  con <- makeContrasts(treat - ctrl, levels=design.matrix)

  #进行差异分析
  qlf <- glmQLFTest(fit, contrast=con)

  dge_results <- as.data.frame(topTags(qlf,n=nrow(dge$counts),sort.by = "logFC"))
  dge_results <- dge_results[order(dge_results$logFC,decreasing = T),]


  dge_results <- dge_results[dge_results$FDR<0.05,]



  ENST_res<- row.names(dge_results)[1:top_n]

  ENST_res


}




M_give_all_AS_factor_logfc <- function(gene_exp,pheno,cluster_to_check){
  AS_list <- read.table("ref/KEGG_SPLICEOSOME.v2024.1.Hs.grp")[-1,]%>%as.vector()



  genelogfc <- M_give_varible_genes_table(gene_exp,pheno,cluster_to_check)

  AS_logfc <- genelogfc[row.names(genelogfc)%in%AS_list,]
  AS_logfc$gene <- row.names(AS_logfc)

  AS.table <- as.data.frame(AS_list)
  colnames(AS.table) <- "gene"

  AS_ALLlogfc <- merge(AS.table,AS_logfc,by="gene",all.x = T)
  AS_ALLlogfc

}



M_plot_visualization_transcript_wholepictures <- function(input_gene_name,
                                                        input_table_dir="output/250306ENST_pct_table_clean.csv",
                                                        pheno="output/250306tsne_kmeans_results.csv",
                                                        group_demo="A",
                                                        anno_list =result,
                                                        color = "#00A087"){
  ##pheno是两列的矩阵，一列是样本名，一列是cluster名
  ##results是列名对应的ENST

  exp <- T_fread_dataframe(input_table_dir)



  pheno <- read_csv(pheno)
  pheno_use <- pheno[pheno$clusters%in%c(group_demo),]
  exp_use <- exp[,colnames(exp)%in%pheno_use$anno.sample]

  ENSG.use <- anno_list$ensembl_gene_id[anno_list$ensembl_transcript_id==input_gene_name]
  ENST.use <- anno_list$ensembl_transcript_id[anno_list$ensembl_gene_id==ENSG.use]
  Symbol.use <- anno_list$external_gene_name[anno_list$ensembl_gene_id==ENSG.use]%>%unique()
  exp_use <- exp_use[row.names(exp_use)%in%ENST.use,]%>%t()%>%as.data.frame()

  #如果一个基因完全没有表达（ENSG上也没有表达的话），其比例无法计算，所以设定为负值，在统计中改为0#
  exp_use[exp_use<0] <- 0

  sapply(exp_use, median)

  exp_use_plot <- data.frame(colnames(exp_use),colMeans(exp_use))
  #exp_use_plot <- data.frame(colnames(exp_use),sapply(exp_use, median))
  colnames(exp_use_plot) <- c("ID","Value")




  ###plot module#############

  p <- ggbarplot(exp_use_plot,
                 x = "ID",
                 y = "Value",
                 fill = color,
                 color = color,
                 palette = NULL,title = paste0(group_demo," ALL ENST in ",ENSG.use," GeneSymbol:",Symbol.use)
  )+ylim(0,1)+coord_flip() # 禁用默认调色板



  p
}









#Tool functions

T_ENST2ENSG<- function(ENST_res){

  ensembl <- readRDS("ref/ensembl_hsp.rds")
  # 获取对应的基因名
  gene_result <- getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"),
                       filters = "ensembl_transcript_id",
                       values = ENST_res,
                       mart = ensembl)

  ENSGres <- gene_result$ensembl_gene_id%>%unique()

}

T_ENST2genename<- function(ENST_res){
  ensembl <- readRDS("ref/ensembl_hsp.rds")
  # 获取对应的基因名
  gene_result <- getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"),
                       filters = "ensembl_transcript_id",
                       values = ENST_res,
                       mart = ensembl)

  ENSGres <- gene_result$external_gene_name%>%unique()




}

T_ENSG2genename<- function(ENSG_res){

  ensembl <- readRDS("ref/ensembl_hsp.rds")
  # 获取对应的基因名
  gene_result <- getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"),
                       filters = "ensembl_gene_id",
                       values = ENSG_res,
                       mart = ensembl)

  ENSGres <- gene_result

}

T_fread_dataframe <- function(dir,col="V1"){
  library(data.table)
  library(magrittr)
  A <- fread(dir)%>%as.data.frame()
  row.names(A)<- A[,col]
  A[,col] <- NULL
  A <- as.data.frame(A)
  A
}


















