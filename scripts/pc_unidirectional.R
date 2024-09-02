required_packages <- c("ggplot2", "tidyr", "dplyr")

# Check which packages are not installed
packages_to_install <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]

# Install the missing packages
if(length(packages_to_install) > 0) {
  install.packages(packages_to_install)
}

# Load the libraries
lapply(required_packages, library, character.only = TRUE)

#library(ggplot2)
#library(tidyr)
#library(dplyr)

folder_path <- "Result/statistics/unidirectional/protein/"

# List all files in the directory
file_list <- list.files(path = folder_path, pattern = "*.tsv", full.names = TRUE)

# Loop through each file and process it
for (file_path in file_list) {

    file_name <- tools::file_path_sans_ext(basename(file_path))


    # Read the TSV file into a data frame
    data2 <- read.table(file_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    #data2
    # Create the input data frame
    my_input <- data.frame(
    my_class = c("gene", "gene", "liftover", "liftover", "Orthologous", "Orthologous", "Orthologous", "One_to_one", "One_to_one", "One_to_one"),
    category = c("pc", "no_pc", "pc_liftover", "no_pc_liftover", 
               "One_to_one_orth", "One_to_many_orth", "One_to_none_orth", 
               "One_pc_to_one_pc", "One_pc_to_one_lncRNA", "One_pc_to_one_other"),
    freq = c(data2$pc, data2$no_pc, data2$pc_liftover, data2$no_pc_liftover, 
           data2$One_to_one_orth, data2$One_to_many_orth, data2$One_to_none_orth, 
           data2$One_pc_to_one_pc, data2$One_pc_to_one_lncRNA, data2$One_pc_to_one_other),
    percentage = c(data2$pc / data2$pc, data2$no_pc_liftover / data2$pc, 
                 data2$pc_liftover / data2$pc, data2$no_pc_liftover / data2$pc, 
                 data2$One_to_one_orth / data2$pc_liftover, data2$One_to_many_orth / data2$pc_liftover, 
                 data2$One_to_none_orth / data2$pc_liftover, 
                 data2$One_pc_to_one_pc / data2$One_to_one_orth, data2$One_pc_to_one_lncRNA / data2$One_to_one_orth, 
                 data2$One_pc_to_one_other / data2$One_to_one_orth),
    conservation = c("conserved", "not_conserved", "conserved", "not_conserved", 
                   "conserved", "conserved", "not_conserved", 
                   "conserved", "conserved", "not_conserved")
    )

    class_levels = c("gene", "liftover", "Orthologous", "One_to_one")
    category_levels = c("pc", "no_pc", 
                    "pc_liftover", "no_pc_liftover", 
                    "One_to_one_orth", "One_to_many_orth", "One_to_none_orth", 
                    "One_pc_to_one_pc", "One_pc_to_one_lncRNA", "One_pc_to_one_other")

    label_levels = c("pc_conserved", "no_pc_not_conserved", 
                 "pc_liftover_conserved", "no_pc_liftover_not_conserved",
                 "One_to_one_orth_conserved", "One_to_many_orth_conserved", "One_to_none_orth_not_conserved",
                 "One_pc_to_one_pc_conserved", "One_pc_to_one_lncRNA_conserved", "One_pc_to_one_other_not_conserved")

    fill_color_vector = c("pc_conserved"="#a15284", 
                      "pc_liftover_conserved"="#a15284",
                      "no_pc_liftover_not_conserved"="#a15284",
                      "One_to_one_orth_conserved"="#a15284",
                      "One_to_many_orth_conserved"="#DD2D4A",
                      "One_to_none_orth_not_conserved"="#74737A",
                      "One_pc_to_one_pc_conserved"="#a15284",
                      "One_pc_to_one_lncRNA_conserved"="#E0C29F",
                      "One_pc_to_one_other_not_conserved"="#ccf8ff")

    alpha_vector = c("pc_conserved"=1,
                 "pc_liftover_conserved"=1,
                 "no_pc_liftover_not_conserved"=0.4,
                 "One_to_one_orth_conserved"=1,
                 "One_to_many_orth_conserved"=0.8,
                 "One_to_none_orth_not_conserved"=0.6,
                 "One_pc_to_one_pc_conserved"=1, 
                 "One_pc_to_one_lncRNA_conserved"=1,
                 "One_pc_to_one_other_not_conserved"=1)

    legend_labels = c("pc_conserved"= paste0("Protein "),
                  "pc_liftover_conserved"= paste0("Lifted"),
                  "no_pc_liftover_not_conserved"= paste0("Non lifted"),
                  "One_to_one_orth_conserved" = paste0("One to one"),
                  "One_to_many_orth_conserved" = paste0("One to many"),
                  "One_to_none_orth_not_conserved"= paste0("One to none"),
                  "One_pc_to_one_pc_conserved" = paste0("One-protein to one-protein"),
                  "One_pc_to_one_lncRNA_conserved" = paste0("One-protein to one-lncRNA"),
                  "One_pc_to_one_other_not_conserved" = paste0("One-protein to one-other"))
    legend_title = paste0("")
    my_input
    # Generate plot input
    my_input$my_class = factor(my_input$my_class, levels=class_levels)
    my_input$category = factor(my_input$category, levels=category_levels)

    # Add label to join the legend
    my_input$label = paste0(my_input$category, "_", my_input$conservation)
    my_input$label = factor(my_input$label, levels=label_levels)

    # Coordinates for the first shaded area
    A_first_x_coord = 1.75
    A_second_x_coord = 1.25
    A_first_y_coord = subset(my_input, my_class=="gene" & category=="no_pc")$percentage

    # Coordinates for the second shaded area
    B_first_x_coord = 2.25
    B_second_x_coord = 2.75
    B_first_y_coord = subset(my_input, my_class=="liftover" & category=="no_pc_liftover")$percentage

    # Coordinates for the third shaded area
    C_first_x_coord = 3.25
    C_second_x_coord = 3.75
    C_first_y_coord = 1 - subset(my_input, my_class=="Orthologous" & category=="One_to_one_orth")$percentage

    # Build dataframe to draw the polygon
    ids = c(rep("A",4), rep("B",4), rep("C",4))
    polygon_dataframe = data.frame(
    id = ids,
    x = c(A_first_x_coord, A_first_x_coord, A_second_x_coord, A_second_x_coord,
        B_first_x_coord, B_first_x_coord, B_second_x_coord, B_second_x_coord,
        C_first_x_coord, C_first_x_coord, C_second_x_coord, C_second_x_coord),
    y = c(A_first_y_coord, 1, 1, 0,
        B_first_y_coord, 1, 1, 0,
        C_first_y_coord, 1, 1, 0)
    )

    # Plot
    compare_plot = ggplot(data=my_input, aes(fill=label)) +
    geom_bar(data=my_input, aes(x=my_class, y=freq, fill=label, alpha=label),
           stat="identity", position="fill", width = 0.5) +
    geom_polygon(data=polygon_dataframe, aes(x=x, y=y, group=id), 
               fill="#156b8a", alpha=0.1) +
    geom_text(data=my_input,
            aes(x=my_class, y=freq, label=ifelse(freq > 0, as.character(freq), "")),
            position=position_fill(vjust=0.5), size=3.5, color="black") +
    scale_fill_manual(values=fill_color_vector, labels=legend_labels) +
    scale_alpha_manual(values=alpha_vector, labels=legend_labels) +
    scale_y_continuous(labels = scales::percent_format(scale = 100)) +
    theme_minimal() +
    theme(axis.title = element_text(color="black", size=10),
        axis.text = element_text(color="black", size=10),
        legend.title = element_blank(),
        legend.position = "right",
        legend.text = element_text(size=8)) +
    labs(x = file_name,
       y = "Percent of genes")

    #print(compare_plot)
    ggsave(paste0("./plots/unidirectional/", file_name, "unidirectional_protein_Orthologous.pdf"), plot = compare_plot, width = 8, height = 4)
}
