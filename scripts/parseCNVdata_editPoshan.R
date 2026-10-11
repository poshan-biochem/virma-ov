# This is a modified parseCNVdata.R script obtained from Garsed group
# Last edited on: 11 Oct 2026

library(data.table)
library(dplyr)
library(plyr)
library(ggplot2)

setwd("../")

cnv_res_file = "data/CnvFacets.cv50_1000.nb500.dp15_5k.ALL.txt.gz"

# read the CNV results
cnv_res = fread(cnv_res_file, sep = "\t", stringsAsFactors = FALSE)

# exclude samples i.e. non-primaries and QC failed
cnv_res <- cnv_res %>% dplyr::filter(! SampleID %in% c("AOCS-076-1-3", "AOCS-166-1-2","MAOC00944-6-X",
                                                       "MAOC01893-6-9","MWMH00552-6-8","MWMH00758-6-3","MWMH00758-8-2"))

# only keep genes of interest
genes_of_interest = c("METTL3","METTL14","METTL16", "METTL5","WTAP", "RBM15", "RBM15B", "ZC3H13", "CBLL1","KIAA1429")
cnv_res <- cnv_res %>% dplyr::filter(GeneID %in% genes_of_interest)

# only keep the segment per sample per gene that is the highest TCN
cnv_res <- cnv_res %>% 
  dplyr::arrange(-`CNV.TCN_EM`) %>% 
  dplyr::group_by(SampleID, GeneID) %>% 
  top_n(n=1, wt = `CNV.TCN_EM`) %>% dplyr::ungroup() %>% data.table()

# create a column that contains a unique segment ID. This will tell us if another gene is in the same segment
cnv_res$SEG_ID <- paste0(paste0("chr",cnv_res$CNV.seqnames),":",cnv_res$CNV.start,"-",cnv_res$CNV.end)

# keep only deletion and amplification events 

df <- cnv_res |> filter(CNV.TCN_EM >=7 | CNV.TCN_EM == 0) |> mutate(CNV.SVTYPE = if_else(CNV.SVTYPE == "DUP-LOH", "DUP", CNV.SVTYPE)) #<- I changed DUP-LOH to DUP. 


# create plot 
plot_df <- as.data.frame(
  table(df$GeneID, df$CNV.SVTYPE))

colnames(plot_df) <- c("GeneID", "CNV.SVTYPE", "n")

plot_df <- plot_df %>%
  filter(n > 0) %>%
  mutate(percentage = n / 126*100) #< 126 is the total number of patients

plot_df$GeneID <- factor(plot_df$GeneID, levels = genes_of_interest)


ggplot(plot_df, aes(
  x = GeneID,
  y = percentage,
  fill = CNV.SVTYPE
)) +
  geom_col(position = "dodge", width = 0.3) +
  labs(
    x = "Gene",
    y = "% of total tumor samples",
    fill = "SV type"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )






