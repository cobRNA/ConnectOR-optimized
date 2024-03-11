library(ggplot2)
library(ggalluvial)
vaccinations <- read.csv('hg38tomouse2.csv')
ggplot(data = vaccinations,
       aes(axis1 = Human, axis2 = Mouse, y = freq)) +
  geom_alluvium(aes(fill = Human),
                curve_type = "arctangent") +
  geom_stratum(color = "grey") +
  #scale_fill_manual(values  = c("darkred", "darkorange4", "darkgoldenrod4", "darkolivegreen4", "cadetblue4")) +
  #scale_color_manual(values = c("darkred", "darkorange4", "darkgoldenrod4", "darkolivegreen4", "cadetblue4")) +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) +
  scale_x_continuous(breaks = 1:2, labels = c("Human", "Mouse"))  +
  theme_minimal() +
    theme(
        legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 14, face = "bold")
    )
