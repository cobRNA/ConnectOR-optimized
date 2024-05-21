library(reshape2)
library(ggplot2)
library(scales)
library(dplyr)

# Reading the data from the CSV file
test2 <- read.table("Result/one_to_one.tsv", header = TRUE, sep = "\t")
melted2 <- melt(test2, id.vars = "Sp1_to_sp2")

melted2$liftover <- ''
melted2 <- melted2 %>%
  mutate(liftover = case_when(
    variable %in% c('pc_to_pc', 'pc_to_other') ~ "protein",
    variable %in% c('lncRNA_to_lncRNA', 'lncRNA_to_other') ~ "lncRNA",
    variable %in% c('ncRNA_to_ncRNA', 'ncRNA_to_other') ~ "ncRNA",

  ))
#'pc_to_pc', 'pc_to_lncRNA', 'pc_to_ncRNA', 'pc_to_other', 'lncRNA_to_lncRNA', 'lncRNA_to_pc', 'lncRNA_to_ncRNA', 'lncRNA_to_other', 'ncRNA_to_ncRNA', 'ncRNA_to_pc', 'ncRNA_to_lncRNA', 'ncRNA_to_other', 'other_to_other', 'other_to_pc', 'other_to_lncRNA', 'other_to_ncRNA'


melted2$liftover <- factor(melted2$liftover, levels = c("protein", "lncRNA", "ncRNA"))
#melted2$Sp1_to_sp2 <- factor(melted2$Sp1_to_sp2, levels = c("hg38tomm39", "mm39tohg38"))
melted2$variable <- factor(melted2$variable, levels = c('pc_to_pc', 'pc_to_other', 'lncRNA_to_lncRNA', 'lncRNA_to_other', 'ncRNA_to_ncRNA', 'ncRNA_to_other'))

p2 <- ggplot(melted2, aes(x = liftover, y = value, fill = variable)) + 
  geom_bar(stat = 'identity', position = 'fill') + 
  geom_text(aes(label = label_comma()(value), y = value), position = position_fill(vjust = 0.5), size = 3) +  
  scale_fill_manual(values = c("pc_to_pc" = "#a15284", 
                                "pc_to_other" = alpha("#a15284", 0.4), 

                                "lncRNA_to_lncRNA" = "#156b8a",
                                "lncRNA_to_other" = alpha("#156b8a", 0.4), 


                                "ncRNA_to_ncRNA" = "#237f5d", 
                                "ncRNA_to_other" = alpha("#237f5d", 0.4)),

                     labels = c("pc_to_pc" = "protein to protein", 
                                "pc_to_other" = "protein to other biotype", 

                                "lncRNA_to_lncRNA" = "lncRNA to lncRNA", 
                                "lncRNA_to_other" = "lncRNA to other biotype", 

                                "ncRNA_to_ncRNA" = "ncRNA to ncRNA",
                                "ncRNA_to_other" = "ncRNA to other biotype")) +   

  facet_wrap(. ~ Sp1_to_sp2, scales="free_y", ncol = 1) +
  labs(x = "Gene biotype", y = "Percent of genes", fill = "Gene's orthology one to one") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(labels = scales::percent_format(scale = 100)) +
  theme_bw()

ggsave("plots/one_to_one.pdf", plot = p2, width = 4.8, height = 6, units = "in", dpi = 600)
ggsave("plots/one_to_one.jpg", plot = p2, width = 4.8, height = 6, units = "in", dpi = 300)

ggsave("plots/one_to_one.png", plot = p2, width = 4.8, height = 6, units = "in", dpi = 900, bg = "transparent")

