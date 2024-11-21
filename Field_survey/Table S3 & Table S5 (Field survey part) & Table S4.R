library(openxlsx)
library(vegan)
# Soil sample grouping information
Field_group = read.xlsx("all_field_data-11-21.xlsx", sheet = "Field_group(all species)11", rowNames = T, colNames = T)
Field_group$Sample_ID = rownames(Field_group)

# Soil sample abundance information
#Field_otu_row2 = read.xlsx("Field_data_asv-11-21.xlsx", sheet = "row_otu", colNames = T, rowNames = T) # row ASVs abundance 
#Field_otu_row = Field_otu_row2[,rownames(Field_group)]
#Field_otu_row[1:6,1:6]
#dim(Field_otu_row)

## checking
#rownames(Field_group) %in% colnames(Field_otu_row)

## Tax INFORMATION
#tax_default <- Field_otu_row2[,(387:394)]

# Resampled by the minimum number of reads per sample
#library(microeco)
#library(mecodev)
#data_default <- microtable$new(sample_table = Field_group,otu_table = Field_otu_row, tax_table = tax_default)
#print(data_default)
#data_default$sample_table
#data_default$tidy_dataset() 
#print(data_default)
#data_default$sample_sums()%>% range

#set.seed(1234)
#data_default$rarefy_samples(sample.size = 9477)
#data_default$sample_sums()%>% range
#print(data_default)
#fungi_Flattening <- data_default$otu_table
#fungi_Flattening[1:5,1:5]

# notes: I have completed the above work, so I directly load the completed file
fungi_Flattening = read.xlsx("fungi_Flattening.xlsx", sheet = "field_flattening", rowNames = T, colNames = T)
fungi_Flattening[1:6, 1:6]
colSums(fungi_Flattening)

################################################################################
########################## Table S3 (Field survey) #############################
################################################################################

Field_group = read.xlsx("all_field_data-11-21.xlsx", sheet = "Field_group(all species)11", rowNames = T, colNames = T)
Field_group$Sample_ID = rownames(Field_group)

# Data Transformation
Field_group$RS = sqrt(Field_group$RS)
Field_group$SRL = log10(Field_group$SRL)
Field_group$Wcont = sqrt(Field_group$Wcont)
Field_group$Soil_N = sqrt(Field_group$Soil_N)
Field_group$Years = as.factor(Field_group$Years)
Field_group$Site = as.factor(Field_group$Site)
Field_group = Field_group[colnames(fungi_Flattening), ] 
unique(Field_group$Species)

# Consider normalizing your data
Field_group_scale = Field_group
shapiro.test(log10(Field_group_scale$Phy_Di))
shapiro.test(log10(Field_group_scale$Fun_Di))
Field_group_scale$Phy_Di_log = log10(Field_group_scale$Phy_Di)
Field_group_scale$Fun_Di_log = log10(Field_group_scale$Fun_Di)

colnames(Field_group_scale)
Field_group_scale[c("Site_pool","Phy_Di","Fun_Di","Phy_Di_log","Fun_Di_log","Soil_ph", "Wcont","Soil_N","Tave","Precipitation","Chol","SLA","LDMC","SRL","FRR","RS")] = 
  scale(Field_group_scale[c("Site_pool","Phy_Di","Fun_Di","Phy_Di_log","Fun_Di_log","Soil_ph", "Wcont","Soil_N","Tave","Precipitation","Chol","SLA","LDMC","SRL","FRR","RS")])

# 
species_hel <- as.data.frame(decostand(t(fungi_Flattening), method = 'hellinger'))
species_hel[1:6,1:6]
dim(species_hel)
BC_dist_field <- vegdist(species_hel, method = 'bray')

rownames(species_hel) %in% rownames(Field_group_scale)

# Table 1 (field survey part)
set.seed(1234)
mod1 = vegan::adonis2(BC_dist_field ~ Years + Site + Family + Species + 
                        Family:Years + Family:Site + Family:Years:Site + 
                        Species:Years + Species:Site,
                      data = Field_group_scale, permutations = 999)
mod1
permanova_TOTAL = as.data.frame(mod1)
permanova_TOTAL$R2 = round(permanova_TOTAL$R2,2)
permanova_TOTAL$`F` = round(permanova_TOTAL$`F`,2)
permanova_TOTAL$`q-vaules`=p.adjust(permanova_TOTAL$`Pr(>F)`, method = "holm")

################################################################################
########################## Table S4 (Field survey) #############################
################################################################################
set.seed(1234)
total_data = Field_group_scale
### full model
library(AICcPermanova)
mod_multifactor99 = with(total_data, adonis2(BC_dist_field ~ Site_pool + Phy_Di_log + Fun_Di_log + Soil_ph + Wcont + Soil_N + Tave + Precipitation + 
                                               Phy_Di_log:Tave + Phy_Di_log:Precipitation + Phy_Di_log:Soil_N + Phy_Di_log:Soil_ph + Phy_Di_log:Wcont + 
                                               Fun_Di_log:Tave + Fun_Di_log:Precipitation + Fun_Di_log:Soil_N + Fun_Di_log:Soil_ph + Fun_Di_log:Wcont, 
                                             data = total_data, permutations = 999, strata = Years))

AICcPermanova::AICc_permanova2(mod_multifactor99)
mod_full_results = as.data.frame(AICcPermanova::AICc_permanova2(mod_multifactor99)); mod_full_results$form = "mod0 ~ full"

## Fun_Di_log:Precipitation
Simplified_mod991 = with(total_data, adonis2(BC_dist_field ~ Site_pool + Phy_Di_log + Fun_Di_log + Soil_ph + Wcont + Soil_N + Tave + Precipitation + 
                                               Phy_Di_log:Tave + Phy_Di_log:Precipitation + Phy_Di_log:Soil_N + Phy_Di_log:Soil_ph + Phy_Di_log:Wcont + 
                                               Fun_Di_log:Tave + Fun_Di_log:Soil_N + Fun_Di_log:Soil_ph + Fun_Di_log:Wcont, 
                                             data = total_data, permutations = 999, strata = Years))

Simplified_mod991
AICcPermanova::AICc_permanova2(Simplified_mod991)
mod_Simplified1 = as.data.frame(AICcPermanova::AICc_permanova2(Simplified_mod991)); mod_Simplified1$form = "mod1 ~ -Fun_Di_log:Precipitation"


## Fun_Di_log:Soil_N
Simplified_mod992 = with(total_data, adonis2(BC_dist_field ~ Site_pool + Phy_Di_log + Fun_Di_log + Soil_ph + Wcont + Soil_N + Tave + Precipitation + 
                                               Phy_Di_log:Tave + Phy_Di_log:Precipitation + Phy_Di_log:Soil_N + Phy_Di_log:Soil_ph + Phy_Di_log:Wcont + 
                                               Fun_Di_log:Tave + Fun_Di_log:Soil_ph + Fun_Di_log:Wcont, 
                                             data = total_data, permutations = 999, strata = Years))

Simplified_mod992
AICcPermanova::AICc_permanova2(Simplified_mod992)
mod_Simplified2 = as.data.frame(AICcPermanova::AICc_permanova2(Simplified_mod992)); mod_Simplified2$form = "mod2 ~ -Fun_Di_log:Soil_N"

## Phy_Di_log:Soil_ph
Simplified_mod993 = with(total_data, adonis2(BC_dist_field ~ Site_pool + Phy_Di_log + Fun_Di_log + Soil_ph + Wcont + Soil_N + Tave + Precipitation + 
                                               Phy_Di_log:Tave + Phy_Di_log:Precipitation + Phy_Di_log:Soil_N + Phy_Di_log:Wcont + 
                                               Fun_Di_log:Tave + Fun_Di_log:Soil_ph + Fun_Di_log:Wcont, 
                                             data = total_data, permutations = 999, strata = Years))

Simplified_mod993
AICcPermanova::AICc_permanova2(Simplified_mod993)
mod_Simplified3 = as.data.frame(AICcPermanova::AICc_permanova2(Simplified_mod993)); mod_Simplified3$form = "mod3 ~ -Phy_Di_log:Soil_ph"


## Fun_Di_log:Tave 
Simplified_mod994 = with(total_data, adonis2(BC_dist_field ~ Site_pool + Phy_Di_log + Fun_Di_log + Soil_ph + Wcont + Soil_N + Tave + Precipitation + 
                                               Phy_Di_log:Tave + Phy_Di_log:Precipitation + Phy_Di_log:Soil_N + Phy_Di_log:Wcont + 
                                               Fun_Di_log:Soil_ph + Fun_Di_log:Wcont, 
                                             data = total_data, permutations = 999, strata = Years))

Simplified_mod994
AICcPermanova::AICc_permanova2(Simplified_mod994)
mod_Simplified4 = as.data.frame(AICcPermanova::AICc_permanova2(Simplified_mod994)); mod_Simplified4$form = "mod4 ~ -Fun_Di_log:Soil_ph"


## Fun_Di_log:Soil_ph 
Simplified_mod995 = with(total_data, adonis2(BC_dist_field ~ Site_pool + Phy_Di_log + Fun_Di_log + Soil_ph + Wcont + Soil_N + Tave + Precipitation + 
                                               Phy_Di_log:Tave + Phy_Di_log:Precipitation + Phy_Di_log:Soil_N + Phy_Di_log:Wcont + 
                                               Fun_Di_log:Wcont, 
                                             data = total_data, permutations = 999, strata = Years))

Simplified_mod995
AICcPermanova::AICc_permanova2(Simplified_mod995)
mod_Simplified5 = as.data.frame(AICcPermanova::AICc_permanova2(Simplified_mod995)); mod_Simplified5$form = "mod5 ~ -Fun_Di_log:Tave"


## Fun_Di_log:Wcont
Simplified_mod996 = with(total_data, adonis2(BC_dist_field ~ Site_pool + Phy_Di_log + Fun_Di_log + Soil_ph + Wcont + Soil_N + Tave + Precipitation + 
                                               Phy_Di_log:Tave + Phy_Di_log:Precipitation + Phy_Di_log:Soil_N + Phy_Di_log:Wcont, 
                                             data = total_data, permutations = 999, strata = Years))

Simplified_mod996
AICcPermanova::AICc_permanova2(Simplified_mod996)
mod_Simplified6 = as.data.frame(AICcPermanova::AICc_permanova2(Simplified_mod996)); mod_Simplified6$form = "mod6 ~ -Phy_Di_log:Soil_N"


## Phy_Di_log:Soil_N
Simplified_mod997 = with(total_data, adonis2(BC_dist_field ~ Site_pool + Phy_Di_log + Fun_Di_log + Soil_ph + Wcont + Soil_N + Tave + Precipitation + 
                                               Phy_Di_log:Tave + Phy_Di_log:Precipitation + Phy_Di_log:Wcont, 
                                             data = total_data, permutations = 999, strata = Years))

Simplified_mod997
AICcPermanova::AICc_permanova2(Simplified_mod997)
mod_Simplified7 = as.data.frame(AICcPermanova::AICc_permanova2(Simplified_mod997)); mod_Simplified7$form = "mod7 ~ -Fun_Di_log:Wcont"

## Phy_Di_log:Precipitation
Simplified_mod998 = with(total_data, adonis2(BC_dist_field ~ Site_pool + Phy_Di_log + Fun_Di_log + Soil_ph + Wcont + Soil_N + Tave + Precipitation + 
                                               Phy_Di_log:Tave + Phy_Di_log:Wcont, 
                                             data = total_data, permutations = 999, strata = Years))

Simplified_mod998
AICcPermanova::AICc_permanova2(Simplified_mod998)
mod_Simplified8 = as.data.frame(AICcPermanova::AICc_permanova2(Simplified_mod998)); mod_Simplified8$form = "mod8 ~ -Phy_Di_log:Precipitation"

### 模型比较
mod_all_results = rbind(mod_full_results,mod_Simplified1,mod_Simplified2,mod_Simplified3,mod_Simplified4,
                        mod_Simplified5,mod_Simplified6,mod_Simplified7,mod_Simplified8)
min_AICc = min(mod_all_results$AICc)
mod_all_results$DeltaAICc = mod_all_results$AICc - min_AICc

### Model Selection-delta_aicc is set to 2
#mod_all_sub = subset(mod_all_results, DeltaAICc < 2)
mod_all_sub = mod_all_results
mod_all_sub$ind_Weight = exp(-0.5 * mod_all_sub$DeltaAICc)
mod_all_sub$AICWeight = mod_all_sub$ind_Weight/sum(mod_all_sub$ind_Weight)
mod_all_sub

################################################################################
# p.adjust
perManova_data = as.data.frame(Simplified_mod998)
perManova_data$`q-vaules`=p.adjust(perManova_data$`Pr(>F)`, method = "holm")
perManova_data$`Pr(>F)` = round(perManova_data$`Pr(>F)`, 3)
perManova_data$`F` = round(perManova_data$`F`, 2)
perManova_data$`q-vaules` = round(perManova_data$`q-vaules`, 3)
perManova_data

################################################################################
########################## Table S5 (Field survey) #############################
################################################################################
# Mantel test
# Single trait
dim(Field_group)
Traits = c("Chol","SLA","LDMC","SRL","FRR","RS")
single_trait_mantel_total = NULL

for (iii in Traits) {
  ## 
  Traits_test = Field_group[,c(iii, "Sample_ID")]
  Traits_test$Sample_ID = NULL
  traits_dis = vegdist(Traits_test, method = 'euclidean')  
  ## BC
  bray_dist = (BC_dist_field)
  set.seed(1234)
  mantel_Bray_Curtis <- vegan::mantel(traits_dis, bray_dist, method = 'spearman', permutations = 999, na.rm = TRUE)
  mantel_Bray_Curtis_result = data.frame(iii,mantel_Bray_Curtis$statistic, mantel_Bray_Curtis$signif)
  colnames(mantel_Bray_Curtis_result) = c("Traits","Mantel_R","p_value")
  single_trait_mantel_total = rbind(single_trait_mantel_total,mantel_Bray_Curtis_result)
}


# all traits data
Traits_test = Field_group[,c("Chol","SLA","LDMC","SRL","FRR","RS","Sample_ID")]
Traits_test$Sample_ID = NULL
traits_dis = vegdist(Traits_test, method = 'euclidean')  

## BC
bray_dist = (BC_dist_field)
set.seed(1234)
mantel_Bray_Curtis <- vegan::mantel(traits_dis, bray_dist, method = 'spearman', permutations = 999, na.rm = TRUE)
all_trait_mantel_total = data.frame(mantel_Bray_Curtis$statistic, mantel_Bray_Curtis$signif)
colnames(all_trait_mantel_total) = c("Mantel_R","p_value")

# Phylogenetic distance
# BC
bray_dist = (BC_dist_field)
##
phylo_dist <- cophenetic(tree)
phylo_dist_data = reshape2::melt(phylo_dist, na.rm = T)
phylo_dist_data$names = paste(phylo_dist_data$Var1, phylo_dist_data$Var2, sep = "_")
colnames(phylo_dist_data)[3] = "Phylo_dist"

# Fungal  distance matrix
fungal_dist_data = reshape2::melt(as.matrix(bray_dist), na.rm = T)
colnames(fungal_dist_data)[1] = "Sample_ID"
fungal_dist_data = fungal_dist_data %>% left_join(Field_group[, c("Species","Sample_ID")], by = "Sample_ID")
colnames(fungal_dist_data)[1] = "Var1"
colnames(fungal_dist_data)[2] = "Sample_ID"
colnames(fungal_dist_data)[3] = "Fungal_dist"
fungal_dist_data = fungal_dist_data %>% left_join(Field_group[, c("Species","Sample_ID")], by = "Sample_ID")
fungal_dist_data$names = paste(fungal_dist_data$Species.x, fungal_dist_data$Species.y, sep = "_")

#
Total_dist_data = fungal_dist_data %>% left_join(phylo_dist_data[, c("names","Phylo_dist")], by = "names")
head(Total_dist_data)
Fungal_dist <- reshape2::dcast(Total_dist_data, Var1  ~ Sample_ID , value.var = "Fungal_dist")
rownames(Fungal_dist) = Fungal_dist$Var1
Fungal_dist = Fungal_dist[,-1]
Phylo_dist <- reshape2::dcast(Total_dist_data, Var1  ~ Sample_ID , value.var = "Phylo_dist")
rownames(Phylo_dist) = Phylo_dist$Var1
Phylo_dist = Phylo_dist[,-1]

colnames(Phylo_dist) %in% colnames(Fungal_dist)
rownames(Phylo_dist) %in% rownames(Fungal_dist)

set.seed(1234)
mantel_Bray_Curtis <- vegan::mantel(as.dist(Phylo_dist), as.dist(Fungal_dist), method = 'spearman', permutations = 999, na.rm = TRUE)
all_phylo_mantel_total = data.frame(mantel_Bray_Curtis$statistic, mantel_Bray_Curtis$signif)
colnames(all_phylo_mantel_total) = c("Mantel_R","p_value")

# single_trait_mantel_total all_phylo_mantel_total all_trait_mantel_total
all_trait_mantel_total$Traits = "Muti_traits"
all_phylo_mantel_total$Traits = "Plant_phylo"
mantel_total_Data = rbind(single_trait_mantel_total, all_trait_mantel_total, all_phylo_mantel_total)
mantel_total_Data$Traits = factor(mantel_total_Data$Traits, levels = c("Plant_phylo","Muti_traits","RS","FRR","SRL","LDMC","SLA","Chol"))
mantel_total_Data$sig <- ifelse(mantel_total_Data$p_value>0.05, '', ifelse(mantel_total_Data$p_value>0.01, '*', ifelse(mantel_total_Data$p_value>0.001, '**', '***')))

# p.adjust
mantel_total_Data$P_value_adj = p.adjust(mantel_total_Data$p_value, method = "holm")
mantel_total_Data$Mantel_R2 <- ifelse(mantel_total_Data$p_value >0.05, '', mantel_total_Data$Mantel_R)
mantel_total_Data$labels = paste(mantel_total_Data$Mantel_R2,"\n",mantel_total_Data$sig, sep = "")
mantel_total_Data$Mantel_R = round(mantel_total_Data$Mantel_R, 2)

library(ggplot2)
mytheme = theme( panel.background = element_rect(fill='white', colour='black'),
                 legend.position = "none",
                 legend.key = element_blank(),
                 #legend.background = element_blank(),   
                 legend.box.background = element_blank(),
                 panel.grid=element_blank(), 
                 legend.title = element_text(size = 9),
                 legend.text = element_text(size = 8),
                 legend.background = element_rect(fill = NA), #axis.ticks.length = unit(0.4,"lines"), 
                 axis.ticks = element_line(color='black'),
                 axis.line = element_line(colour = "black"), 
                 axis.title.x = element_text(colour='black', size=13),
                 axis.title.y = element_text(colour='black', size=13),
                 axis.text = element_text(colour='black',size=11),
                 plot.tag = element_text(size = 14, face = "bold")) 

ggplot(mantel_total_Data, aes(Traits,"BC_dist")) +
  geom_tile(aes(fill=Mantel_R), color = "white") +
  geom_text(aes(label=labels), color="black", size=3) + 
  scale_fill_gradient2(low = 'white',mid = 'white',high ='#f99104')+
  labs(x=NULL,y=NULL,fill = NULL, tag = "d") + 
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
                              "Muti_traits" = "Trait dissmilarity", "Plant_phylo" = "Phylogenetic distance"),position = "top")



### Notice that,
### For more picture details, we have further adjusted it in Adobe illustrator.



