library(reshape2)
library(ggplot2)
library(scales)
library(dplyr)

# Reading the data from the CSV file
test2 <- read.table("Result/otholog.tsv", header = TRUE, sep = "\t")
melted2 <- melt(test2, id.vars = "Sp1_to_sp2")
    #print(total_pc, one_to_one_pc, one_to_many_pc,one_to_none_pc)                  
    #print(total_lncRNA, one_to_one_lncRNA, one_to_many_lncRNA,one_to_none_lncRNA)
    #print(total_ncRNA, one_to_one_ncRNA, one_to_many_ncRNA,one_to_none_ncRNA)
    #print(total_other, one_to_one_other, one_to_many_other,one_to_none_other)
# Assigning liftover values based on conditions using dplyr
melted2 <- melted2 %>%
  mutate(liftover = case_when(
    variable %in% c('one_to_one_pc', 'one_to_many_pc', 'one_to_none_pc') ~ "protein",
    variable %in% c('one_to_one_lncRNA', 'one_to_many_lncRNA', 'one_to_none_lncRNA') ~ "lncRNA",
    variable %in% c('one_to_one_ncRNA', 'one_to_many_ncRNA', 'one_to_none_ncRNA') ~ "ncRNA",
  ))

# Reorder liftover factor levels
melted2$liftover <- factor(melted2$liftover, levels = c("protein", "lncRNA", "ncRNA"))

# Combine intergenic categories and update other categories
melted2 <- melted2 %>%
  mutate(variable = case_when(
    variable %in% c('one_to_none_pc', 'one_to_none_lncRNA', 'one_to_none_ncRNA') ~ "one_to_none",
    variable %in% c('one_to_many_pc', 'one_to_many_lncRNA', 'one_to_many_ncRNA') ~ "Intergenic",
    variable %in% c('one_to_one_pc') ~ "one_to_one_pc",
    variable %in% c('one_to_one_lncRNA') ~ "one_to_one_lncRNA",
    variable %in% c('one_to_one_ncRNA') ~ "one_to_one_ncRNA"

  ))

# Reorder Sp1_to_sp2 factor levels
#melted2$Sp1_to_sp2 <- factor(melted2$Sp1_to_sp2, levels = c("hg38tomm39", "mm39tohg38"))
melted2$variable <- factor(melted2$variable, levels = c('one_to_one_pc', 'one_to_one_lncRNA', 'one_to_one_ncRNA', 'Intergenic', 'one_to_none'))
# Define colors and labels
color_vector <- c("one_to_one_pc" = "#a15284", 
                  "one_to_one_lncRNA" = "#156b8a", 
                  "one_to_one_ncRNA" = "#237f5d", 
                  "Intergenic" = "#f0be39", 
                  "one_to_none" = "#7f8c8d")

label_vector <- c("one_to_one_pc" = "one to one", 
                  "one_to_one_lncRNA" = "one to one ", 
                  "one_to_one_ncRNA" = "one to one ",
                  "Intergenic" = "one to many",
                  "one_to_none" = "one to none")

# Update the label with a thousands separator
p2 <- ggplot(melted2, aes(x = liftover, y = value, fill = variable)) + 
  geom_bar(stat = 'identity', position = 'fill') + 
  geom_text(aes(label = comma(value)), position = position_fill(vjust = 0.5), size = 3) + 
  scale_fill_manual(values = color_vector, labels = label_vector) + 
  facet_wrap(. ~ Sp1_to_sp2, scales="free_y", ncol = 1) +
  labs(x = "Gene biotype", y = "Percent of genes", fill = "Gene's orthology type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(labels = scales::percent_format(scale = 100)) +
  theme_bw() 

# Saving the plot as PDF and TIFF
ggsave("plots/ortholog.pdf", plot = p2, width = 4.6, height = 6, units = "in", dpi = 600)
ggsave("plots/ortholog.jpg", plot = p2, width = 4.6, height = 6, units = "in", dpi = 300)
ggsave("plots/ortholog.png", plot = p2, width = 4.6, height = 6, units = "in", dpi = 900, bg = "transparent")

