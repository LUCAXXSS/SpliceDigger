
SpliceDigger.findMarkers <- function(SD_object=A,use_data="iso_exp_pct_filtered",
                                     anno_to_use="sort",cluster_to_check="LCBM",
                                     only.pos=T,logfc.threshold=0.5,FDR.threshold=0.05
                                     ){

  exp_deg <- SD_object[[use_data]]




  group = SD_object$anno[[anno_to_use]]


  group[!group==cluster_to_check] <- "ctrl"
  group[group==cluster_to_check] <- "treat"

  group <- factor(group,levels = c("ctrl","treat"))

  dge <- edgeR::DGEList(counts = exp_deg, group = group)
  # 重新计算库大小
  dge <- edgeR::calcNormFactors(dge)

  # 估计共同离散度和趋势离散度
  dge <- edgeR::estimateDisp(dge)

  design.matrix<-model.matrix(~0+group)
  colnames(design.matrix)<-levels(group)

  fit <- edgeR::glmQLFit(y = dge, design = design.matrix)

  con <- limma::makeContrasts(treat - ctrl, levels=design.matrix)

  #进行差异分析
  qlf <- edgeR::glmQLFTest(fit, contrast=con)

  dge_results <- as.data.frame(edgeR::topTags(qlf,n=nrow(dge$counts),sort.by = "logFC"))
  dge_results <- dge_results[order(dge_results$logFC,decreasing = T),]

  if (only.pos) {
    dge_results <- dge_results[dge_results$logFC>logfc.threshold&dge_results$FDR<FDR.threshold,]
  }
  dge_results

}
