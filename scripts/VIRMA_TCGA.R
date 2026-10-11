# TCGA analysis of OC samples

library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(ggplot2)

query <- GDCquery(
  project = "TCGA-OV",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(query, directory = "TCGA-OV", files.per.chunk = 10)

ov <- GDCprepare(query, directory = "TCGA-OV")

saveRDS(ov, file = "TCGA-OV/TCGA_OV")

ov <- readRDS("~/TCGA-OV/TCGA_OV")

ov_coldata <- colData(ov) |> as_tibble()

# plot VIRMA expression versus CNA events ----

virma_expression <- read_tsv("~/Downloads/plot.txt")

virma_expression <- janitor::clean_names(virma_expression)

virma_expression$copy_number_alterations |> table()

virma_expression$copy_number_alterations <- factor(virma_expression$copy_number_alterations, 
                                                   levels = c("Not Profiled",
                                                              "VIRMA: Shallow Deletion",
                                                              "VIRMA: Diploid",
                                                              "VIRMA: Gain",
                                                              "VIRMA: Amplification"))

ggplot(virma_expression, aes(x = copy_number_alterations, y = virma_m_rna_expression_rna_seq_v2_rsem_log2_value_1_3)) + geom_boxplot()

