################################################################################
########################### Fig. S1 (Field survey) #############################
################################################################################
library(openxlsx)
library(nlme)
library(car)
library(MuMIn)
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

Field_group = read.xlsx("all_field_data-11-21.xlsx", sheet = "Field_group(all species)11", rowNames = T, colNames = T)
Field_group$Sample_ID = rownames(Field_group)
#Field_group$Years = as.factor(Field_group$Years)
Field_group$Site = factor(Field_group$Site, levels = unique(Field_group$Site[order(Field_group$Latitude)]))

#### Environmental effects
Field_group$Effect_size = log(Field_group$Fungal_field_Di/Field_group$Fungal_green_Di)

# Data Transformation
Field_group$RS = sqrt(Field_group$RS)
Field_group$SRL = log10(Field_group$SRL)
Field_group$Wcont = sqrt(Field_group$Wcont)
Field_group$Soil_N = sqrt(Field_group$Soil_N)
Field_group$Phy_Di_log = log10(Field_group$Phy_Di)
Field_group$Fun_Di_log = log10(Field_group$Fun_Di)
unique(Field_group$Species)

# Consider normalizing your data
pd_attributes_variable <- attributes(scale(Field_group[c("Site_pool","Phy_Di","Fun_Di","Phy_Di_log","Fun_Di_log","Soil_ph", "Wcont","Soil_N","Tave","Precipitation")]))

total_data = Field_group
colnames(total_data)
total_data[c("Site_pool","Phy_Di","Fun_Di","Phy_Di_log","Fun_Di_log","Soil_ph", "Wcont","Soil_N","Tave","Precipitation")] = 
  scale(total_data[c("Site_pool","Phy_Di","Fun_Di","Phy_Di_log","Fun_Di_log","Soil_ph", "Wcont","Soil_N","Tave","Precipitation")])

fm1 = lme(Fungal_SR ~ Site_pool + Phy_Di_log + Fun_Di_log + Soil_ph + Wcont + Soil_N + Precipitation + Tave +
            Phy_Di_log:Tave + Phy_Di_log:Precipitation + Phy_Di_log:Soil_ph + Phy_Di_log:Wcont + Phy_Di_log:Soil_N + 
            Fun_Di_log:Tave + Fun_Di_log:Precipitation + Fun_Di_log:Soil_ph + Fun_Di_log:Wcont + Fun_Di_log:Soil_N, 
          random = ~1|Years, method = "ML", data = total_data)

# Simplify the model - use Log-likelihood ration test backed up by AIC
drop1(fm1, test = "Chi") # 
f.sbs1 <- update(fm1, ~. -Fun_Di_log:Wcont)
AIC(f.sbs1)
drop1(f.sbs1, test = "Chi")
f.sbs2 <- update(f.sbs1, ~. -Phy_Di_log:Soil_ph)
AIC(f.sbs2)
drop1(f.sbs2, test = "Chi")
f.sbs3 <- update(f.sbs2, ~. -Phy_Di_log:Wcont)
AIC(f.sbs3)
drop1(f.sbs3, test = "Chi")
f.sbs4 <- update(f.sbs3, ~. -Phy_Di_log:Precipitation)
AIC(f.sbs4)
drop1(f.sbs4, test = "Chi")
f.sbs5 <- update(f.sbs4, ~. -Phy_Di_log:Soil_N)
AIC(f.sbs5)
drop1(f.sbs5, test = "Chi")
f.sbs6 <- update(f.sbs5, ~. -Wcont)
AIC(f.sbs6)
drop1(f.sbs6, test = "Chi")
f.sbs7 <- update(f.sbs6, ~. -Fun_Di_log:Soil_N)
AIC(f.sbs7)
drop1(f.sbs7, test = "Chi")
f.sbs8 <- update(f.sbs7, ~. -Fun_Di_log:Tave)
drop1(f.sbs8, test = "Chi")
f.sbs9 <- update(f.sbs8, ~. -Phy_Di_log:Tave)
drop1(f.sbs9, test = "Chi")
f.sbs10 <- update(f.sbs9, ~. -Tave)
drop1(f.sbs10, test = "Chi")
f.sbs11 <- update(f.sbs10, ~. -Phy_Di_log )
drop1(f.sbs11, test = "Chi")
f.sbs12 <- update(f.sbs11, ~. -Fun_Di_log:Precipitation)
drop1(f.sbs12, test = "Chi")
f.sbs13 <- update(f.sbs12, ~. -Fun_Di_log:Soil_ph)
AIC(f.sbs13)
drop1(f.sbs13, test = "Chi")
f.sbs14 <- update(f.sbs13, ~. -Soil_ph)
AIC(f.sbs14)
drop1(f.sbs14, test = "Chi")

library(lmtest)
lrtest(f.sbs4,f.sbs5,f.sbs6,f.sbs7,f.sbs8,f.sbs9,f.sbs10,f.sbs11,f.sbs12,f.sbs13,f.sbs14)
anova(f.sbs4,f.sbs5,f.sbs6,f.sbs7,f.sbs8,f.sbs9,f.sbs10,f.sbs11,f.sbs12,f.sbs13,f.sbs14)
AIC_order = as.data.frame(AIC(f.sbs4,f.sbs5,f.sbs6,f.sbs7,f.sbs8,f.sbs9,f.sbs10,f.sbs11,f.sbs12,f.sbs13,f.sbs14));AIC_order

f.sbs.fin <- update(f.sbs14, ~., method = "REML") # Update to REML to extract estimates
shapiro.test(residuals(f.sbs.fin))
summary(f.sbs.fin)
Anova(f.sbs.fin)
MuMIn::r.squaredGLMM(f.sbs.fin)
as.data.frame(vif(f.sbs.fin))
summary(f.sbs.fin)

### p.adjust
tables2<- as.data.frame(car::Anova(f.sbs.fin, type = 3))
tables2=tables2[-which(rownames(tables2) == "(Intercept)"),]
tables2=as.data.frame(tables2)
tables2$`q-vaules`=p.adjust(tables2$`Pr(>Chisq)`, method = "holm")
tables2$`p-value` = round(tables2$`Pr(>Chisq)`, 3)
tables2$`q-vaules` = round(tables2$`q-vaules`, 3)
tables2$Parameter = rownames(tables2)

library(effectsize)
MegaModelSummary = as.data.frame(effectsize(f.sbs.fin))[-1,]
MegaModelSummary$Parameter2 = c("Site pool","Funct-Dist","Soil N", "Precipitation")
MegaModelSummary$Parameter2 = factor(MegaModelSummary$Parameter2, levels = rev(c("Site pool","Funct-Dist","Soil N", "Precipitation")))
MegaModelSummary = MegaModelSummary %>% left_join(tables2[,c("Parameter","q-vaules")], by = "Parameter")
MegaModelSummary$p_label = paste0("p = ",round(MegaModelSummary$`q-vaules`, 3))
MegaModelSummary$sig2 <- as.factor(ifelse(MegaModelSummary$`q-vaules` > 0.05, 2, 1))

MegaModelSummary$Group = c("Site pool", "Plant attributes", "Soil properties", "Climate")
MegaModelSummary$Group = factor(MegaModelSummary$Group, levels = c("Site pool", "Plant attributes", "Soil properties", "Climate", "Interaction"))

ggplot(MegaModelSummary, aes(x = Parameter2, y = Std_Coefficient, fill = Group))+
  geom_errorbar(aes(ymin=CI_low, ymax=CI_high), width=0.1, size = 0.8, color = "black")+
  geom_point(size = 3.5, pch = 21)+
  geom_segment(aes(y = 0, yend = 0, x = 0.5, xend = 4.35), color = "black", linetype = "dashed") + 
  geom_text(aes(y = CI_high, label = p_label),            
            hjust =  -0.4, vjust = 0.2, angle = 0, color = "black",size = 3.2) + 
  #annotate("text", x = 4.5 , y = 0.08,label = "Best model: R2m = 0.08, R2c = 0.12", colour="black", size = 4) +  
  annotate("text", x = 4.5, y = 0.08,
           label = expression("Best model: " * italic(R)^2 * "m = 0.08, " * italic(R)^2 * "c = 0.12"),
           colour = "black", size = 4) + 
  labs(x = '', y = 'Standard regression coefficients', color = '') +  
  theme_classic() + coord_flip() +  
  scale_fill_manual(values = c("Site pool" = "#949698","Plant attributes" = "#DE7963", "Soil properties" = "#3F425A",
                               "Climate" = "#82B19D","Interaction" = "#F0C986")) +
  theme(axis.text.x = element_text(color = "black", size = 11),
        axis.text.y = element_text(color = "black", size = 11),
        axis.title.x =  element_text(color = "black", size = 13),
        legend.text = element_text(size = 9, color = "black"),
        plot.margin = margin(0.5,1.5,0.5,1.5, unit = "cm"),
        legend.position = 'none',legend.title = element_blank(),legend.key = element_blank(),
        plot.tag = element_text(size = 14, face = "bold")) +
  scale_shape_manual(values = c(16,21)) -> p_a; p_a


#
library(glmm.hp)
hierarchical_data = as.data.frame(glmm.hp(f.sbs.fin, type = "R2")$hierarchical.partitioning)
hierarchical_data$Group = c("Site pool", "Plant attributes", "Soil properties", "Climate")
hierarchical_data$Parameter2 = c("Site pool","Funct-Dist", "Soil N", "Precipitation")
hierarchical_data$Group = factor(hierarchical_data$Group, levels = c("Site pool", "Plant attributes", "Soil properties", "Climate", "Interaction"))

hierarchical_data2 = hierarchical_data %>% group_by(Group) %>% 
  summarise(`I.perc(%)` = sum(`I.perc(%)`))

ggplot()+
  geom_bar(data = hierarchical_data2, aes(x = "", y = `I.perc(%)`, fill = Group), 
           stat = "identity", width = 0.5, color = "black")+
  scale_y_continuous(expand = c(0, 0), position = "right")+
  scale_fill_manual(values = c("Site pool" = "#949698","Plant attributes" = "#DE7963", "Soil properties" = "#3F425A",
                               "Climate" = "#82B19D","Interaction" = "#F0C986")) +
  theme_classic()+
  theme(legend.background = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.text.y = element_text(color = "black", size = 11),
        axis.title.y = element_text(color = "black", size = 13),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_line(linewidth = NA))+
  guides(fill = guide_legend(title = "Province",ncol = 3, nrow = 2, override.aes = list(shape = 22, size=0.2))) +
  labs(y = "Relative effect of estimates (%)") -> p_b; p_b

library(patchwork)
(p_a+p_b) + plot_layout(widths = c(0.7,0.3)) -> P_s1; P_s1
