library(tidyverse)
library(janitor)

df <- read_tsv("~/Downloads/CNA_Genes.txt") # CNA_Genes.txt was downloaded from cBIoPortal

df <- janitor::clean_names(df)
df$freq <- gsub("%", replacement = "",df$freq)
df$freq <- as.numeric(df$freq)

genes <- c("METTL3", "METTL14", "WTAP", "VIRMA", "YTHDF1", "YTHDF2", "YTHDF3", "YTHDC1", "YTHDC2", "FTO", "ALKBH5")

df_mtc <- df |> filter(df$gene %in% genes) |> select(gene, freq, cna)


ggplot(df_mtc, aes(x = gene, y = freq, fill = cna)) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.4
  ) +
  labs(
    x = NULL,
    y = "% of samples with alterations",
    fill = "Condition"
  ) +
  theme_classic(base_size = 14)


df_mtc <- df_mtc[order(df_mtc$gene),]


write.csv(df_mtc, file = "mtc_cna_analysis.csv")


# Figure -2 density plot

library(dplyr)
library(purrr)
library(ggplot2)

# read files 
files <- list.files("~/Downloads", "*OC.csv", full.names = T)

crispr_effects <- lapply(files, FUN = 'read.csv')

files

names(crispr_effects) <- c("ALKBH5", "FTO", "METTL14", "METTL3", "VIRMA", "WTAP", "YTHDC1", "YTHDC2", "YTHDF1", "YTHDF2", "YTHDF3")

names(crispr_effects)

# convert data into longer format for ease of plotting 

crispr_long <- map2_df(
  crispr_effects,
  names(crispr_effects),
  ~{
    score_col <- grep("Chronos", names(.x), value = TRUE)
    
    .x %>%
      transmute(
        Gene = .y,
        Score = .data[[score_col]],
        Disease = OncotreePrimaryDisease,
        Lineage = OncotreeLineage,
        Cell_id = `Cell.Line.Name`
      )
  }
)

# add factors to change the order later
crispr_long$Gene <- factor(
  crispr_long$Gene,
  levels = c("METTL3", "METTL14", "WTAP", "VIRMA", "YTHDF1", "YTHDF2", "YTHDF3", "YTHDC1", "YTHDC2", "FTO", "ALKBH5")
)

head(crispr_long)

# calculate the mean lethality score 
means <- crispr_long %>% group_by(Gene) %>% 
summarise(mean_score = mean(Score, na.rm = TRUE))

# generate a density plot
ggplot(crispr_long,
       aes(x = Score,
           fill = Gene,
           colour = Gene)) +
  geom_density(alpha = 0.2) +
  geom_rug(
    sides = "b",
    alpha = 1,
    length = unit(0.08, "npc")
  ) +
  geom_vline(
    data = means,
    aes(xintercept = mean_score),
    colour = "red",
    linetype = "dashed",
    linewidth = 0.8,
    inherit.aes = FALSE
  ) +
  facet_grid(Gene ~ .) +
  theme_classic()

ggsave("crispr_effects.tiff", dpi = 600, height = 12, width = 5)



