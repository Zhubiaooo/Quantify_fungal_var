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

# Loading R packages
library(openxlsx)
library(vegan)
library(dplyr)
library(phytools)
library(funrar)

# Read in the fungal abundance information table
Green_otu <- read.xlsx("Greenhouse_data_asv.xlsx", sheet = "Overall_otu", colNames = T, rowNames = T)
Green_otu[1:6,1:6]
Green_otu <- Green_otu[,-c(1)]
dim(Green_otu)

# Hellinger Transformation
species_hel <- as.data.frame(decostand(t(Green_otu), method = 'hellinger'))
species_hel[1:6,1:6]

# Read in sample grouping information
Green_group <- read.xlsx("Greenhouse_data_asv.xlsx", sheet = "all_otu_data2", colNames = T, rowNames = T)
Green_group$sample <- rownames(Green_group)
colnames(Green_group)

# Read plant functional traits (mean values)
traits_mean <- read.xlsx("traits_mean.xlsx", sheet = "traits_mean", colNames = T, rowNames = T)
colnames(traits_mean)
# Transformation
shapiro.test(log10(traits_mean$RS))
shapiro.test(log10(traits_mean$SRL))
traits_mean$RS = sqrt(traits_mean$RS)
traits_mean$SRL = log10(traits_mean$SRL)

# Permutational multivariate analysis of variance to explore the effects of plant families and species on fungal composition
set.seed(1234)
bray_dist <- vegdist(species_hel, method = 'bray')
adonis_result <- adonis2(bray_dist ~ Family + Species, Green_group, permutations = 999)
adonis_result
# p.adjust
permanova_data2 = as.data.frame(adonis_result)
permanova_data2$`q-vaules`=p.adjust(permanova_data2$`Pr(>F)`, method = "holm")

### Visualization (PCoA)
pcoa <- cmdscale(bray_dist , k = 2, eig = TRUE)
plot_data <- data.frame({pcoa$points})[1:2]
plot_data$sample <- rownames(plot_data)
names(plot_data)[1:2] <- c('PCoA1', 'PCoA2')
eig = pcoa$eig
plot_data <- merge(plot_data, Green_group, by = 'sample', all.x = TRUE)
head(plot_data)[1:6]

plot_data$Family = as.factor(plot_data$Family)
df = plot_data %>% group_by(Family) %>% 
  summarise(PCoA1_mean = mean(PCoA1), PCoA1_se = sd(PCoA1)/(sqrt(length(PCoA1))),
            PCoA2_mean = mean(PCoA2), PCoA2_se = sd(PCoA2)/(sqrt(length(PCoA2))))
#Rmisc::summarySE(plot_data, measurevar = c("PCoA1"), groupvars = c("Years", "Latitude", "Origin"))

plot_data2 = merge(plot_data, df, by = c("Family"))
unique(df$Family)
disss = 0.08
ggplot(df, aes(PCoA1_mean, PCoA2_mean))+
  geom_point(plot_data2,mapping = aes(PCoA1, PCoA2, color = Family), size = 2, pch= 16) +
  geom_segment(plot_data2,mapping = aes(x = PCoA1_mean, y = PCoA2_mean, xend = PCoA1, yend = PCoA2, color = Family),
               linetype = 1, linewidth = 0.2, alpha = 0.4)+
  geom_errorbar(data = df,mapping = aes(ymax = PCoA2_mean+PCoA2_se, ymin=PCoA2_mean-PCoA2_se, color = Family),width=0.01,size=0.3,alpha = 1)+#
  geom_errorbarh(data = df,mapping = aes(xmax=PCoA1_mean+PCoA1_se,xmin=PCoA1_mean-PCoA1_se, color = Family),height=0.01,size=0.3,alpha = 1) +
  geom_point(df,mapping = aes(PCoA1_mean, PCoA2_mean, color = Family, fill = Family), size = 4, pch = 21, color = "black")+
  scale_color_manual(values = family_color) +
  scale_fill_manual(values = family_color) +
  scale_shape_manual(values = c(21,16)) +
  labs(x=paste("PCoA1 (", format(100 * eig[1] / sum(eig), digits=3), "%)", sep=""),
       y=paste("PCoA2 (", format(100 * eig[2] / sum(eig), digits=3), "%)", sep=""),
       tag = "a")+
  #geom_vline(aes(xintercept = 0),linetype="dashed")+
  #geom_hline(aes(yintercept = 0),linetype="dashed")+
  annotate("text", x = -0.05, y = -0.16, label = paste0("Family: R2 =", round(adonis_result$R2[1],2),", p = ", adonis_result$`Pr(>F)`[1])) + 
  annotate("text", x = -0.05, y = -0.20, label = paste0("Species: R2 =", round(adonis_result$R2[2],2),", p = ", adonis_result$`Pr(>F)`[2])) + 
  Pcoa_theme +
  theme(plot.tag = element_text(size = 14, face = "bold")) + 
  scale_x_continuous(labels = scales::label_comma(accuracy =0.01))+
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01)) -> Pcoa_all

library(ggmagnify)
library(ggfx)
ggm <- ggmagnify(Pcoa_all,
                 xlim = c(-0.15, -0.03), ylim = c(-0.08, 0.08),
                 inset_xlim = c(-0.20, 0.18), inset_ylim = c(0.24, 0.45)); ggm

# Relationship between single traits and trait matrix, phylogenetic distance matrix and community composition (Mantel test)
trait_names = c("Chol","SLA","LDMC","SRL","FRR","RS")

# single traits
single_trait_mantel = NULL
Green_group2 = Green_group %>% left_join(traits_mean, by = c("Species", "Origin"))
colnames(Green_group2)

cor.test(Green_group2$Overall_Richness, Green_group2$SRL)
cor.test(Green_group2$Overall_Richness, Green_group2$LDMC)

for (i in trait_names) {
  all_otu_TEST = Green_group2[,c(i, "Species")]
  all_otu_TEST$Species = NULL
  traits_dis <- vegdist(all_otu_TEST, method = 'euclidean')  
  ##Bray-Curtis
  set.seed(1234)
  mantel_Bray_Curtis <- vegan::mantel(traits_dis, bray_dist, method = 'spearman', permutations = 999, na.rm = TRUE)
  mantel_Bray_Curtis_result = data.frame("Test",mantel_Bray_Curtis$statistic, mantel_Bray_Curtis$signif,i,'Bray_Curtis')
  colnames(mantel_Bray_Curtis_result) = c("Test","Mantel_R","P_value","Traits","dist_type")
  single_trait_mantel = rbind(single_trait_mantel,mantel_Bray_Curtis_result)
}
single_trait_mantel

# Trait matrix
Fun_dist <- compute_dist_matrix(traits_mean[,c("Chol","SLA","LDMC","FRR","SRL","RS")], metric = "euclidean", scale = TRUE, center = TRUE)
Fun_dist_data = reshape2::melt(Fun_dist, na.rm = T)
Fun_dist_data$names = paste(Fun_dist_data$Var1, Fun_dist_data$Var2, sep = "_")
colnames(Fun_dist_data)[3] = "Funct_dist"
head(Fun_dist_data)

# Phylogenetic distance matrix
tree <- read.newick("iq_tree.NEWICK")
plot(tree)
to_drop<-c("Amborella_trichopoda","")
tree <- drop.tip(as.phylo(tree), to_drop) 
plant_dist <- cophenetic(tree)
plant_dist_data = reshape2::melt(plant_dist, na.rm = T)
plant_dist_data$names = paste(plant_dist_data$Var1, plant_dist_data$Var2, sep = "_")
colnames(plant_dist_data)[3] = "Phylo_dist"
head(plant_dist_data)

Root_dist_data = reshape2::melt(as.matrix(bray_dist), na.rm = T)
colnames(Root_dist_data)[1] = "sample"
Root_dist_data = Root_dist_data %>% left_join(all_otu[, c("Species","sample")], by = "sample")
colnames(Root_dist_data)[1] = "Var1"
colnames(Root_dist_data)[2] = "sample"
colnames(Root_dist_data)[3] = "Fungal_dist"
Root_dist_data = Root_dist_data %>% left_join(all_otu[, c("Species","sample")], by = "sample")
Root_dist_data$names = paste(Root_dist_data$Species.x, Root_dist_data$Species.y, sep = "_")
head(Root_dist_data)

#
Total_dist_data = Root_dist_data %>% left_join(plant_dist_data[, c("names","Phylo_dist")], by = "names") %>% 
  left_join(Fun_dist_data[, c("names","Funct_dist")], by = "names")

head(Total_dist_data)

Fungal_dist <- reshape2::dcast(Total_dist_data, Var1  ~ sample , value.var = "Fungal_dist")
rownames(Fungal_dist) = Fungal_dist$Var1
Fungal_dist = Fungal_dist[,-1]
##
Phylo_dist <- reshape2::dcast(Total_dist_data, Var1  ~ sample , value.var = "Phylo_dist")
rownames(Phylo_dist) = Phylo_dist$Var1
Phylo_dist = Phylo_dist[,-1]
##
Funct_dist <- reshape2::dcast(Total_dist_data, Var1  ~ sample , value.var = "Funct_dist")
rownames(Funct_dist) = Funct_dist$Var1
Funct_dist = Funct_dist[,-1]

colnames(Phylo_dist) %in% colnames(Fungal_dist)
rownames(Phylo_dist) %in% rownames(Fungal_dist)

rownames(Funct_dist) %in% rownames(Fungal_dist)
rownames(Funct_dist) %in% rownames(Fungal_dist)

set.seed(1234)
mantel_fun <- vegan::mantel(as.dist(Funct_dist), as.dist(Fungal_dist), method = 'spearman', permutations = 999, na.rm = TRUE)
mantel_result_fun = data.frame("Test",mantel_fun$statistic, mantel_fun$signif,"Muti_traits",'Bray_Curtis')
colnames(mantel_result_fun) = c("Test","Mantel_R","P_value","Traits","dist_type")

set.seed(1234)
mantel_phy <- vegan::mantel(as.dist(Phylo_dist), as.dist(Fungal_dist), method = 'spearman', permutations = 999, na.rm = TRUE)
mantel_result_phy = data.frame("Test",mantel_phy$statistic, mantel_phy$signif,"Plant_phylo",'Bray_Curtis')
colnames(mantel_result_phy) = c("Test","Mantel_R","P_value","Traits","dist_type")

## 
mantel_total_Data = rbind(single_trait_mantel, mantel_result_fun, mantel_result_phy)
mantel_total_Data$Traits = factor(mantel_total_Data$Traits, levels = c("Plant_phylo","Muti_traits","RS","FRR","SRL","LDMC","SLA","Chol"))
mantel_total_Data$sig <- ifelse(mantel_total_Data$P_value>0.05, '', ifelse(mantel_total_Data$P_value>0.01, '*', ifelse(mantel_total_Data$P_value>0.001, '**', '***')))

select_traits = c("Chol","SLA","LDMC","FRR","SRL","RS","Muti_traits","Plant_phylo")
mantel_total_Data$Mantel_R = round(mantel_total_Data$Mantel_R, 2)

mantel_total_Data$Mantel_R2 <- ifelse(mantel_total_Data$P_value >0.05, '', mantel_total_Data$Mantel_R)
mantel_total_Data$labels = paste(mantel_total_Data$Mantel_R2,"\n",mantel_total_Data$sig, sep = "")

col <- colorRampPalette(c("#00d1b2","white","#f99104"))(50)

p = ggplot(mantel_total_Data, aes(Traits,dist_type)) +
  geom_tile(aes(fill=Mantel_R), color = "white") +
  #geom_tile(fill = "white",color = "white")+
  #geom_point(aes(fill = Mantel_R),color = "white",shape = 22, size = 20)+ 
  geom_text(aes(label=labels), color="black", size=3) + 
  scale_fill_gradient2(low = 'white',mid = 'white',high ='#f99104')+
  #scale_fill_gradient2(low = "#5C5DAF",mid = "white",high = "#EA2E2D") +
  labs(x=NULL,y=NULL,fill = NULL, tag = "b") + 
  mytheme + 
  theme(panel.grid = element_blank(), 
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background = element_rect(color = 'black', fill = 'transparent'),
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size=11,color = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "bottom") +
  coord_flip() +
  guides(fill = guide_colorbar(barwidth = 8, barheight = 0.8,
                               title.theme = element_text(size = 12),
                               label.theme = element_text(size = 10),
                               frame.colour = "black", frame.linewidth = 0.01)) +
  scale_x_discrete(labels = c("Chol" = "Leaf chlorophyll", "SLA" = "Specific leaf area",
                              "LDMC" = "Leaf dry matter content","SRL" = "Specific root length",
                              "FRR" = "Fine-root mass ratio", "RS" = "Root-shoot ratio",
                              "Muti_traits" = "Euclidean trait distance", "Plant_phylo" = "Plant phylogenetic distance"), position = "top")

p

library(grid)
p + annotation_custom(
  grob = linesGrob(gp = gpar(col = "black", lwd = 1)), 
  ymin = 0.0, ymax = 8, xmin = 2.5, xmax = 2.5) -> pp

# p.adjust
mantel_total_Data$P_value_adj = p.adjust(mantel_total_Data$P_value, method = "holm")
