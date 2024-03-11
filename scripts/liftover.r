library(reshape2)
library(ggplot2)
library(scales)
library(dplyr)


test2 <- read.table("Result/liftover.tsv", header = TRUE, sep = "\t")
melted2 <- melt(test2, id.vars = "Sp1_to_sp2")
melted2$liftover <- ''
melted2 <- melted2 %>%
  mutate(liftover = case_when(
    variable %in% c('protein_lifted', 'non_lifted_protein') ~ "protein",
    variable %in% c('lncRNA_lifted', 'non_lifted_lncRNA') ~ "lncRNA",
    variable %in% c('ncRNA_lifted', 'non_lifted_ncRNA') ~ "ncRNA",
    TRUE ~ NA_character_  # Closing the case_when function properly
  ))
melted2$liftover <- factor(melted2$liftover, levels = c("protein", "lncRNA", "ncRNA"))
melted2$Sp1_to_sp2 <- factor(melted2$Sp1_to_sp2, levels = c("hg38tomm39", "hg38todanrer11", "hg38todm6", "mm39tohg38", "mm39todanrer11", "mm39todm6", "danrer11tohg38", "danrer11tomm39", "danrer11todm6", "dm6todanrer11", "dm6tomm39", "dm6tohg38"))
melted2$variable <- factor(melted2$variable, levels = c('non_lifted_protein', 'protein_lifted', 'non_lifted_lncRNA', 'lncRNA_lifted', 'non_lifted_ncRNA', 'ncRNA_lifted'))
p2 <- ggplot(melted2, aes(x = liftover, y = value, fill = variable)) + 
  geom_bar(stat = 'identity', position = 'fill') + 
  geom_text(aes(label = label_comma()(value), y = value), position = position_fill(vjust = 0.5), size = 3) +  
  scale_fill_manual(values = c("protein_lifted" = "#a15284", 
                                "non_lifted_protein" = alpha("#a15284", 0.2), 
                                "lncRNA_lifted" = "#156b8a", 
                                "non_lifted_lncRNA" = alpha("#156b8a", 0.2), 
                                "ncRNA_lifted" = "#237f5d", 
                                "non_lifted_ncRNA" = alpha("#237f5d", 0.2)),
                     labels = c("protein_lifted" = "protein lifted", 
                                "non_lifted_protein" = "non lifted protein", 
                                "lncRNA_lifted" = "lncRNA lifted", 
                                "non_lifted_lncRNA" = "non lifted lncRNA", 
                                "ncRNA_lifted" = "ncRNA lifted", 
                                "non_lifted_ncRNA" = "non lifted ncRNA")) +   
  facet_wrap(. ~ Sp1_to_sp2, scales="free_y", ncol = 3) +
  labs(x = "Gene biotype", y = "Percent of genes", fill = "Gene's orthology type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(labels = scales::percent_format(scale = 100)) +
  theme_bw()

ggsave("plots/lifted_gene.pdf", plot = p2, width = 8, height = 8, units = "in", dpi = 600)




