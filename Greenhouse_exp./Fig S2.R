# Custom Style
library(ggplot2)
mytheme = theme( panel.background = element_rect(fill='white', colour='black'),
                 legend.position = "none",
                 panel.grid=element_blank(), 
                 legend.title = element_blank(),
                 legend.text = element_text(size = 8),
                 legend.background = element_rect(fill = NA), #axis.ticks.length = unit(0.4,"lines"), 
                 axis.ticks = element_line(color='black'),
                 axis.line = element_line(colour = "black"), 
                 axis.title.x = element_text(colour='black', size=13),
                 axis.title.y = element_text(colour='black', size=13),
                 axis.text = element_text(colour='black',size=11),
                 plot.tag = element_text(size = 14, face = "bold")) 
family_color = c("Acanthaceae"="#332500","Amaranthaceae"="#542D20","Asteraceae" ="#994240","Caryophyllaceae" = "#D75C4D", "Cyperaceae" = "#E68C51",
                 "Euphorbiaceae" ="#F59D52", "Fabaceae" = "#EFBA55","Lamiaceae" ="#FCD170","Malvaceae"="#FEE1B0", "Onagraceae" = "#C5E0D0",
                 "Phytolaccaceae" = "#ABDCE0","Poaceae" ="#7FC5D8","Polygonaceae"="#73BCD5", "Solanaceae" = "#528FAC", "Urticaceae" = "#376694", "Verbenaceae" = "#1F466F")

library(openxlsx)
all_otu <- read.xlsx("Greenhouse_data_asv.xlsx", sheet = "all_otu_data2", colNames = T, rowNames = T)

#### Effect of family and species on fungal richness
colnames(all_otu)
mod = lm(Overall_Richness ~ Family + Species , data = all_otu)
tables2 = as.data.frame(anova(mod))
tables2$`q-vaules`=p.adjust(tables2$`Pr(>F)`, method = "holm")
tables2$`Pr(>F)` = round(tables2$`Pr(>F)`, 4)
tables2$`q-vaules` = round(tables2$`q-vaules`, 4)

species_order <- read.xlsx("Greenhouse_data_asv.xlsx", sheet = "species_list", colNames = T, rowNames = F)
colnames(species_order) = c("order","Species2","Species")
all_otu$Species2 <- gsub("_", " ", all_otu$Species)
all_otu$Species2 <- factor(all_otu$Species2, levels = species_order$Species2)
all_otu$sig = ifelse(all_otu$Origin == "Native"," ","*")
all_otu$Species3 = paste0(all_otu$Species2,all_otu$sig)
colnames(all_otu)
all_otu$Sample_ID = rownames(all_otu)
all_otu2 = all_otu %>% left_join(species_order[,c(1,3)], by = "Species")

# Fig S2a
ggplot(all_otu2,aes(x = Family, y = Overall_Richness, fill = Family)) + 
  geom_boxplot(width = 0.6,alpha = 0.75) + 
  #scale_color_manual(values = family_color) +
  scale_fill_manual(values = family_color) +
  mytheme + 
  theme(axis.text.y = element_text(face = "italic",size=10)) + 
  coord_flip() + 
  labs(x = NULL, y = "Overall fungal richness", tag = "a") -> p1; p1

# Fig S2b
ggplot(all_otu2,aes(x = reorder(Species3, -order), y = Overall_Richness, fill = Family)) + 
  geom_boxplot(width = 0.6,alpha = 0.75) + 
  #scale_color_manual(values = family_color) +
  scale_fill_manual(values = family_color) +
  mytheme + 
  theme(axis.text.y = element_text(face = "italic",size=10)) + 
  coord_flip() + 
  labs(x = NULL, y = "Overall fungal richness", tag = "b") -> p2; p2

library(patchwork)
p1|p2





