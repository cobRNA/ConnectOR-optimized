library(reshape2)
library(ggplot2)
library(scales)
library(dplyr)
#library(forcats)
# Reading the data from the CSV file
test2 <- read.table("Result/ortholog.tsv", header = TRUE, sep = "\t")
melted2 <- melt(test2, id.vars = "Sp1_to_sp2")

# Assigning liftover values based on conditions using dplyr
melted2 <- melted2 %>%
  mutate(liftover = case_when(
    variable %in% c('protein_ortholog', 'protein_lncRNA_ncRNA_ortholog', 'protein_intergenic') ~ "protein",
    variable %in% c('lncRNA_ortholog', 'lncRNA_protein_ncRNA_ortholog', 'lncRNA_intergenic') ~ "lncRNA",
    variable %in% c('ncRNA_ortholog', 'ncRNA_protein_lncRNA_ortholog', 'ncRNA_intergenic') ~ "ncRNA"
  ))

# Reorder liftover factor levels
melted2$liftover <- factor(melted2$liftover, levels = c("protein", "lncRNA", "ncRNA"))

# Reorder Sp1_to_sp2 factor levels
melted2$Sp1_to_sp2 <- factor(melted2$Sp1_to_sp2, levels = c("hg38tomm39", "hg38todanrer11", "hg38todm6", "mm39tohg38", "mm39todanrer11", "mm39todm6", "danrer11tohg38", "danrer11tomm39", "danrer11todm6", "dm6todanrer11", "dm6tomm39", "dm6tohg38"))

# Update the label with a thousands separator
p2 <- ggplot(melted2, aes(x = liftover, y = value, fill = variable)) + 
  geom_bar(stat = 'identity', position = 'fill') + 
  geom_text(aes(label = comma(value)), position = position_fill(vjust = 0.5), size = 3.2) +  # Removed fill aesthetic

  scale_fill_manual(values = c("protein_ortholog" = "#a15284", 
                                "protein_lncRNA_ncRNA_ortholog" = "#FF9671", 
                                "protein_intergenic" = "#bbcfd8", 
                                "lncRNA_ortholog" = "#156b8a", 
                                "lncRNA_protein_ncRNA_ortholog" = "#FF9671", 
                                "lncRNA_intergenic" = "#bbcfd8", 
                                "ncRNA_ortholog" = "#237f5d", 
                                "ncRNA_protein_lncRNA_ortholog" = "#FF9671", 
                                "ncRNA_intergenic" = "#bbcfd8"),
                   labels = c("protein_ortholog" = "Protein Ortholog", 
                              "protein_intergenic" = "Protein Intergenic", 
                              "lncRNA_ortholog" = "lncRNA Ortholog", 
                              "lncRNA_protein_ncRNA_ortholog" = "lncRNA ortholog with other biotypes", 
                              "lncRNA_intergenic" = "lncRNA Intergenic", 
                              "ncRNA_ortholog" = "ncRNA Ortholog", 
                              "ncRNA_protein_lncRNA_ortholog" = "ncRNA ortholog with other biotypes", 
                              "ncRNA_intergenic" = "ncRNA Intergenic")) + 
  facet_wrap(. ~ Sp1_to_sp2, scales="free_y", ncol = 3) +
  labs(x = "Gene biotype", y = "Percent of genes", fill = "Gene's orthology type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(labels = scales::percent_format(scale = 100)) +
  theme_bw()

# Saving the plot as PDF and TIFF
ggsave("plots/ortholog.pdf", plot = p2, width = 10, height = 8, units = "in", dpi = 600)



