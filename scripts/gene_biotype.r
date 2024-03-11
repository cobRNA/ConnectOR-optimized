library(reshape2)
library(ggplot2)
library(scales)

# Read data
#test2 <- read.csv("gene_stat.csv", header = TRUE, stringsAsFactors = FALSE)
test2 <- read.table("Result/gene_stat.tsv", header = TRUE, sep = "\t")
# Melt the data
melted2 <- melt(test2, "Species")
melted2$Species <- factor(melted2$Species, levels = c("hg38", "mm39", "danrer11", "dm6"))
# Plot with facet_wrap
p <- ggplot(melted2, aes(x = variable, y = value, fill = variable)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = comma(value), y = value), color = "black", size = 3, vjust = 0.99, position = position_dodge(width = 0.9)) +  
  scale_fill_manual(values = c("protein" = "#a15284", "lncRNA" = "#156b8a", "ncRNA" = "#237f5d", "other" = "#CADBC0")) +
  facet_wrap(~ Species, scales = "free_y", ncol = 2) +
  labs(x = "Gene biotype", y = "Number of Genes") +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_bw() +  # Set background to white
  guides(fill = "none")  

# Save the plot
ggsave("plots/gene_biotype_plot.pdf", plot = p, width = 8, height = 8, units = "in", dpi = 600)



