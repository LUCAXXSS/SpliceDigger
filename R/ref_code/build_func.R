




rm(list = ls())
library(biomaRt)

# 连接到 Ensembl 数据库
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# 查询 transcript_id (ENST) 和对应的 gene_id (ENSG)
enst_to_ensg <- getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"),
                      mart = ensembl)

enst_ensg_genename <- enst_to_ensg


save(enst_ensg_genename,file = "data/enst_ensg_genename.rda")

load("data/enst_ensg_genename.rda")

# 保存到本地文件
write.csv(enst_to_ensg, "ENST_to_ENSG_mapping.csv", row.names = FALSE)


A <- fread("R/ref_code/EXP_notsoundTPM_250306.csv")%>%as.data.frame()
row.names(A)<- A[,"V1"]
A[,"V1"] <- NULL
A <- as.data.frame(A)
A

exp_You_etal_2025 <- A


usethis::use_data(exp_You_etal_2025, overwrite = TRUE)

ENST <- row.names(iso_tpm_You_etal_2025)

data_ta <- .T_ENST2TABLE(ENST)


devtools::load_all()

A <- SpliceDigger.create_SD_object(ENST_table = iso_tpm_You_etal_2025,exp_table= exp_You_etal_2025)


anno <- read.csv("~/ymy/240524AS_LCBM_analysis/input/250214pheno_detailed.csv",row.names = 1)%>%as.data.frame()

anno <- anno[row.names(anno)%in%colnames(exp_You_etal_2025),,drop = FALSE]


anno <- anno[match(colnames(exp_You_etal_2025),row.names(anno)),,drop = FALSE]

dataanno_You_etal_2025 <- anno

#usethis::use_data(dataanno_You_etal_2025)

A <- SpliceDigger.add_annotation(A,dataanno_You_etal_2025)

A <- SpliceDigger.filter_SDobject(AS_threshold = 0.2,AS_var_topn = 5000,exp_threshold = 20,exp_var_topn = 5000)

A <- SpliceDigger.calculate_tSNE(A,use_data = "iso_exp_pct_filtered",tSNE_perplexity = 15,kmeans_centers = 3,seed = 4401)

SpliceDigger.plot_tSNE(plot_dataset ="iso_exp_pct_filtered",anno_to_show = "iso_exp_pct_filtered_kmeans_cluster")


df <- SpliceDigger.findMarkers(anno_to_use = "iso_exp_pct_filtered_kmeans_cluster",cluster_to_check = "cluster_3",only.pos = T)


SpliceDigger.visualize_AS_transcript(SD_object = A,anno_to_use = "iso_exp_pct_filtered_kmeans_cluster",
                                     anno_to_demo = "cluster_3",feature = "SEC61G",
                                     feature_type = "Symbol",
                                     calculate_method = "mean",color = "#DC0000")

SpliceDigger.visualize_exp_heatmap(SD_object = A,use_data = "exp_raw",FDR_threshold = 0.05,
                                   anno_to_use = "sort",feature_to_plot = c("EGFR","PDGFRA","SRSF1","SRSF2","SRSF3","SRSF4"))



