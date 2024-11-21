################################################################################
############################ Fig 3 (Field survey) ##############################
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

fm1 = lme(Effect_size ~ Site_pool + Phy_Di_log + Fun_Di_log + Soil_ph + Wcont + Soil_N + Precipitation + Tave +
            Phy_Di_log:Tave + Phy_Di_log:Precipitation + Phy_Di_log:Soil_ph + Phy_Di_log:Wcont + Phy_Di_log:Soil_N + 
            Fun_Di_log:Tave + Fun_Di_log:Precipitation + Fun_Di_log:Soil_ph + Fun_Di_log:Wcont + Fun_Di_log:Soil_N, 
          random = ~1|Years, method = "ML", data = total_data)

# Simplify the model - use Log-likelihood ration test backed up by AIC
drop1(fm1, test = "Chi") # 
f.sbs1 <- update(fm1, ~. -Phy_Di_log:Soil_N)
AIC(f.sbs1)
drop1(f.sbs1, test = "Chi")
f.sbs2 <- update(f.sbs1, ~. -Phy_Di_log:Soil_ph)
AIC(f.sbs2)
drop1(f.sbs2, test = "Chi")
f.sbs3 <- update(f.sbs2, ~. -Phy_Di_log:Wcont)
AIC(f.sbs3)
drop1(f.sbs3, test = "Chi")
f.sbs4 <- update(f.sbs3, ~. -Soil_ph)
AIC(f.sbs4)
drop1(f.sbs4, test = "Chi")
f.sbs5 <- update(f.sbs4, ~. -Fun_Di_log:Soil_ph)
AIC(f.sbs5)
drop1(f.sbs5, test = "Chi")
f.sbs6 <- update(f.sbs5, ~. -Phy_Di_log:Precipitation)
AIC(f.sbs6)
drop1(f.sbs6, test = "Chi")
f.sbs7 <- update(f.sbs6, ~. -Fun_Di_log:Wcont)
AIC(f.sbs7)
drop1(f.sbs7, test = "Chi")
f.sbs8 <- update(f.sbs7, ~. -Fun_Di_log:Precipitation)
drop1(f.sbs8, test = "Chi")
f.sbs9 <- update(f.sbs8, ~. -Precipitation)
drop1(f.sbs9, test = "Chi")
f.sbs10 <- update(f.sbs9, ~. -Wcont)
drop1(f.sbs10, test = "Chi")

library(lmtest)
lrtest(f.sbs4,f.sbs5,f.sbs6,f.sbs7,f.sbs8,f.sbs9,f.sbs10)
anova(f.sbs4,f.sbs5,f.sbs6,f.sbs7,f.sbs8,f.sbs9,f.sbs10)
AIC_order = as.data.frame(AIC(f.sbs4,f.sbs5,f.sbs6,f.sbs7,f.sbs8,f.sbs9,f.sbs10));AIC_order
round(MuMIn::Weights(AIC(f.sbs4,f.sbs5,f.sbs6,f.sbs7,f.sbs8,f.sbs9,f.sbs10)), 3)


f.sbs.fin <- update(f.sbs9, ~., method = "REML") # Update to REML to extract estimates

#f.sbs.fin = lme(Effect_size ~ Site_pool + Phy_Di_log + Fun_Di_log + Wcont + Soil_N + Tave + 
#                  Phy_Di_log:Tave + Fun_Di_log:Tave + Fun_Di_log:Soil_N, random = ~1|Years, data = total_data)

shapiro.test(residuals(f.sbs.fin))
summary(f.sbs.fin)
car::Anova(f.sbs.fin)
MuMIn::r.squaredGLMM(f.sbs.fin)
as.data.frame(vif(f.sbs.fin))


### p.adjust
tables2<-summary(f.sbs.fin)$tTable
tables2=tables2[-which(rownames(tables2) == "(Intercept)"),]
tables2=as.data.frame(tables2)
tables2$`q-vaules`=p.adjust(tables2$`p-value`, method = "holm")
tables2$`p-value` = round(tables2$`p-value`, 3)
tables2$`q-vaules` = round(tables2$`q-vaules`, 3)
tables2$Parameter = rownames(tables2)


library(flextable)
library(dplyr)
# needs the packages flextable and officer to be installed

my_theme_flextable=function (x) 
{
  if (!inherits(x, "flextable")) 
    stop("my_theme supports only flextable objects.")
  x <- border_remove(x)
  x <- hline_top(x, part = "header", border = officer::fp_border())
  x <- hline_bottom(x, part = "header", border = officer::fp_border())
  x <- hline_bottom(x, part = "body", border = officer::fp_border())
  x <- font(x,fontname = "Times")
  x <- font(x,part="header",fontname = "Times")
  x
}


table_anova=data.frame(Step=c("0","1","2","3","4","5","6","7","8","9","10"),
                       Model=c("Full model","-Phylo Di × Soil pH","-Phylo Di × Soil N","-Phylo Di × Wcont","-Funct Di × Soil pH","-Soil pH","-Phylo Di × Precipitation","-Funct Di × Temperature","-Funct Di × Wcont", "-Wcont", "-Precipitation"),
                       #n_par=c(anova(fm1,f.sbs1)[1,1],anova(fm1,f.sbs1)[2,1],anova(f.sbs1,f.sbs2)[2,1],anova(f.sbs2,f.sbs3)[2,1],anova(f.sbs3,f.sbs4)[2,1],anova(f.sbs4,f.sbs5)[2,1],anova(f.sbs5,f.sbs6)[2,1],anova(f.sbs6,f.sbs7)[2,1],anova(f.sbs7,f.sbs8)[2,1],anova(f.sbs8,f.sbs9)[2,1],anova(f.sbs9,f.sbs10)[2,1]),
                       Chisq=round(c(NA,anova(fm1,f.sbs1)[2,8],anova(f.sbs1,f.sbs2)[2,8],anova(f.sbs2,f.sbs3)[2,8],anova(f.sbs3,f.sbs4)[2,8],anova(f.sbs4,f.sbs5)[2,8],anova(f.sbs5,f.sbs6)[2,8],anova(f.sbs6,f.sbs7)[2,8],anova(f.sbs7,f.sbs8)[2,8],anova(f.sbs8,f.sbs9)[2,8],anova(f.sbs9,f.sbs10)[2,8]),1),
                       df=round(c(NA,anova(fm1,f.sbs1)[2,3],anova(f.sbs1,f.sbs2)[2,3],anova(f.sbs2,f.sbs3)[2,3],anova(f.sbs3,f.sbs4)[2,3],anova(f.sbs4,f.sbs5)[2,3],anova(f.sbs5,f.sbs6)[2,3],anova(f.sbs6,f.sbs7)[2,3],anova(f.sbs7,f.sbs8)[2,3],anova(f.sbs8,f.sbs9)[2,3],anova(f.sbs9,f.sbs10)[2,3]),1),
                       p_value=as.character(signif(c(NA,anova(fm1,f.sbs1)[2,9],anova(f.sbs1,f.sbs2)[2,9],anova(f.sbs2,f.sbs3)[2,9],anova(f.sbs3,f.sbs4)[2,9],anova(f.sbs4,f.sbs5)[2,9],anova(f.sbs5,f.sbs6)[2,9],anova(f.sbs6,f.sbs7)[2,9],anova(f.sbs7,f.sbs8)[2,9],anova(f.sbs8,f.sbs9)[2,9],anova(f.sbs9,f.sbs10)[2,9]),3)),
                       R2m=round(c(r.squaredGLMM(fm1)[1],r.squaredGLMM(f.sbs1)[1],r.squaredGLMM(f.sbs2)[1],r.squaredGLMM(f.sbs3)[1],r.squaredGLMM(f.sbs4)[1],r.squaredGLMM(f.sbs5)[1],r.squaredGLMM(f.sbs6)[1],r.squaredGLMM(f.sbs7)[1],r.squaredGLMM(f.sbs8)[1],r.squaredGLMM(f.sbs9)[1],r.squaredGLMM(f.sbs10)[1]),2),
                       R2c=round(c(r.squaredGLMM(fm1)[2],r.squaredGLMM(f.sbs1)[2],r.squaredGLMM(f.sbs2)[2],r.squaredGLMM(f.sbs3)[2],r.squaredGLMM(f.sbs4)[2],r.squaredGLMM(f.sbs5)[2],r.squaredGLMM(f.sbs6)[2],r.squaredGLMM(f.sbs7)[2],r.squaredGLMM(f.sbs8)[2],r.squaredGLMM(f.sbs9)[2],r.squaredGLMM(f.sbs10)[2]),2),
                       AIC=round(c(AIC(fm1),AIC(f.sbs1),AIC(f.sbs2),AIC(f.sbs3),AIC(f.sbs4),AIC(f.sbs5),AIC(f.sbs6),AIC(f.sbs7),AIC(f.sbs8),AIC(f.sbs9),AIC(f.sbs10)),1),
                       AIC_weight=round(MuMIn::Weights(AIC(fm1,f.sbs1,f.sbs2,f.sbs3,f.sbs4,f.sbs5,f.sbs6,f.sbs7,f.sbs8,f.sbs9,f.sbs10)), 3)) %>%
  mutate(delta_AIC=round(AIC-min(AIC),1)) %>%
  dplyr::select(Step:AIC,delta_AIC,AIC_weight) %>%
  flextable() %>%
  bold(i=~delta_AIC==0) %>%
  set_header_labels(Chisq="χ²",p_value="p-value",R2c="R²c",R2m="R²m",delta_AIC="∆AIC",AIC_weight="AIC weight" ) %>%
  my_theme_flextable() %>%
  autofit()
table_anova

#
library(effectsize)
MegaModelSummary = as.data.frame(effectsize(f.sbs.fin))[-1,]
MegaModelSummary$Parameter2 = c("Site pool","Phylo-Dist", "Funct-Dist", "Wcont", "Soil N", "Temperature", 
                                "Phylo-Dist × Temperature", "Funct-Dist × Temperature", "Funct-Dist × Soil N")
MegaModelSummary$Parameter2 = factor(MegaModelSummary$Parameter2, levels = rev(c("Site pool","Phylo-Dist", "Funct-Dist", "Wcont", "Soil N", "Temperature", 
                                                                                 "Phylo-Dist × Temperature", "Funct-Dist × Temperature", "Funct-Dist × Soil N")))
MegaModelSummary = MegaModelSummary %>% left_join(tables2[,c("Parameter","q-vaules")], by = "Parameter")
MegaModelSummary$p_label = paste0("p = ",round(MegaModelSummary$`q-vaules`, 3))
MegaModelSummary$sig2 <- as.factor(ifelse(MegaModelSummary$`q-vaules` > 0.05, 2, 1))

MegaModelSummary$Group = c("Site pool", "Plant attributes", "Plant attributes", "Soil properties", "Soil properties",
                           "Climate", "Interaction", "Interaction", "Interaction")
MegaModelSummary$Group = factor(MegaModelSummary$Group, levels = c("Site pool", "Plant attributes", "Soil properties", "Climate", "Interaction"))



ggplot(MegaModelSummary, aes(x = Parameter2, y = Std_Coefficient, fill = Group))+
  geom_errorbar(aes(ymin=CI_low, ymax=CI_high), width=0.2, size = 0.8, color = "black")+
  geom_point(size = 3.5, pch = 21)+
  geom_segment(aes(y = 0, yend = 0, x = 0.5, xend = 9.2), color = "black", linetype = "dashed") + 
  geom_text(aes(y = CI_high, label = p_label),            
            hjust =  -0.4, vjust = 0.2, angle = 0, color = "black",size = 3.2) + 
  #annotate("text", x = 9.5 , y = 0.08,label = "Best model: R2m = 0.30, R2c = 0.32", colour="black", size = 4) +  
  annotate("text", x = 9.46, y = 0.08,
           label = expression("Best model: " * italic(R)^2 * "m = 0.30, " * italic(R)^2 * "c = 0.32"),
           colour = "black", size = 4) + 
  labs(x = '', y = 'Standard regression coefficients', color = '', tag = "a") +  
  theme_classic() + coord_flip() +  
  scale_fill_manual(values = c("#949698","#DE7963", "#3F425A","#82B19D","#F0C986")) +
  theme(axis.text.x = element_text(color = "black", size = 11),
        axis.text.y = element_text(color = "black", size = 11),
        axis.title.x =  element_text(color = "black", size = 13),
        legend.text = element_text(size = 9, color = "black"),
        plot.margin = margin(0.5,1.5,0.5,1.5, unit = "cm"),
        legend.position = 'none',legend.title = element_blank(),legend.key = element_blank(),
        plot.tag = element_text(size = 14, face = "bold")) +
  scale_shape_manual(values = c(16,21)) -> p3a; p3a

library(glmm.hp)
hierarchical_data = as.data.frame(glmm.hp(f.sbs.fin, type = "R2")$hierarchical.partitioning)
hierarchical_data$Group = c("Site pool", "Plant attributes", "Plant attributes", "Soil properties", "Soil properties",
                            "Climate", "Interaction", "Interaction", "Interaction")
hierarchical_data$Parameter2 = c("Site pool","Phylo-Dist", "Funct-Dist", "Wcont", "Soil N", "Temperature", 
                                 "Phylo-Dist × Temperature", "Funct-Dist × Temperature", "Funct-Dist × Soil N")
hierarchical_data$Group = factor(hierarchical_data$Group, levels = c("Site pool", "Plant attributes", "Soil properties", "Climate", "Interaction"))

hierarchical_data2 = hierarchical_data %>% group_by(Group) %>% 
  summarise(`I.perc(%)` = sum(`I.perc(%)`))

ggplot()+
  geom_bar(data = hierarchical_data2, aes(x = "", y = `I.perc(%)`, fill = Group), 
           stat = "identity", width = 0.5, color = "black")+
  scale_y_continuous(expand = c(0, 0), position = "right")+
  scale_fill_manual(values = c("#949698","#DE7963", "#3F425A","#82B19D","#F0C986")) +
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
  labs(y = "Relative effect of estimates (%)") -> p3b; p3b


# Phylo-Dist × Temperature
library(effects)
eff_mod<-effect("Phy_Di_log:Tave",f.sbs.fin,xlevels=5)
eff_mod_data <- data.frame(eff_mod)

# backtransform
eff_mod_data["Tave"] <- pd_attributes_variable$`scaled:center`["Tave"] + pd_attributes_variable$`scaled:scale`["Tave"]*eff_mod_data["Tave"]
eff_mod_data["Phy_Di_log"] <- pd_attributes_variable$`scaled:center`["Phy_Di_log"] + pd_attributes_variable$`scaled:scale`["Phy_Di_log"]*eff_mod_data["Phy_Di_log"]

total_data["Phy_Di_row"] <- pd_attributes_variable$`scaled:center`["Phy_Di_log"] + pd_attributes_variable$`scaled:scale`["Phy_Di_log"]*total_data["Phy_Di_log"]
total_data["Tave_row"] <- pd_attributes_variable$`scaled:center`["Tave"] + pd_attributes_variable$`scaled:scale`["Tave"]*total_data["Tave"]
eff_mod_data$Phy_Di_log = round(eff_mod_data$Phy_Di_log, 2)
ggplot()+
  #geom_point(data = total_data, mapping = aes(x = Tave_row, y = Effect_size), pch = 21) + 
  geom_line(data = eff_mod_data, mapping = aes(Tave,fit,color=factor(Phy_Di_log)), size=1.25)+
  labs(x = expression("Annual temperature (°C)"), 
       y = "Effect size", tag = "b", fill = expression("Phylo Di(log"[10]*"(10)")) +
  geom_ribbon(data = eff_mod_data, mapping = aes(x = Tave, ymin=fit-se,ymax=fit+se,fill=factor(Phy_Di_log)),alpha=0.3, colour= NA)+
  #scale_fill_manual(values=c("#FEEC00","#C0B667", "#737891", "#425C8E", "#001B56"),guide=F)+
  #scale_color_manual(values=c("#FEEC00","#C0B667", "#737891", "#425C8E", "#001B56"),breaks=waiver(), name="Phylo Di (log10)")+
  scale_fill_manual(values=c("#184C3F","#6CB2AA", "#E4CB8F", "#BF8E43", "#57320F"),guide=F)+
  scale_color_manual(values=c("#184C3F","#6CB2AA", "#E4CB8F", "#BF8E43", "#57320F"),breaks=waiver(), name="Phylo Di (log10)")+
  guides(color=guide_legend(override.aes=list(fill=NA)))+
  annotate("text", label = expression(italic(p) == 0.024), x = 15, y = 0.25, size = 4) + 
  mytheme  + 
  theme(legend.position = c(0.8,0.28))-> P3b; P3b

library(patchwork)
(p3a+p3b) + plot_layout(widths = c(0.7,0.3)) -> P3a

(P3a|P3b) + plot_layout(widths = c(0.6,0.4))

### Notice that,
### For more picture details, we have further adjusted it in Adobe illustrator.
