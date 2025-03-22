rm(list = ls())




library(tidyverse)
library(data.table)

# 加载 biomaRt
library(biomaRt)

########导入基因列表############

exp <- fread("input/iso_tpm_formatted_1.txt")





# 选择人类基因数据库
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# 定义 ENST 列表
enst_ids <- exp$V1  # 替换为你的 ENST ID 列表

# 获取对应的基因名
result <- getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"),
                filters = "ensembl_transcript_id",
                values = enst_ids,
                mart = ensembl)

#对每一个标本都计算一次

exp_merged <- merge(result,exp,by.x = "ensembl_transcript_id",by.y = "V1")

exp_merged_dict <- exp_merged[,colnames(exp_merged)%in%c("ensembl_transcript_id","ensembl_gene_id")]

exp_merged <- exp_merged[,!colnames(exp_merged)%in%c("ensembl_transcript_id","external_gene_name")]

#########把每一个基因都计算好重复的比例，by gpt###########
exp_merged_result <- exp_merged %>%
  group_by(ensembl_gene_id) %>%
  mutate(across(where(is.numeric), ~ . / sum(.), .names = "pct_{.col}"))

exp_merged_result <- exp_merged_result[, grepl("pct", names(exp_merged_result))]%>%as.data.frame()

row.names(exp_merged_result) <- exp_merged_dict$ensembl_transcript_id

write.csv(exp_merged_result,file = "output/ENST_pct_table_raw.csv",row.names = T)

#把所有NaN替换成-1值 
exp_merged_result_clean <- exp_merged_result %>%
  mutate(across(everything(), ~ ifelse(is.nan(.), -1, .)))

exp_merged_result_clean <- exp_merged_result_clean[!rowMeans(exp_merged_result_clean)==0,]
exp_merged_result_clean <- exp_merged_result_clean[!rowMeans(exp_merged_result_clean)==1,]

write.csv(exp_merged_result_clean,file = "output/ENST_pct_table_clean.csv",row.names = T)
