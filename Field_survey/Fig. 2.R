################################################################################
############################ Fig 2 (Field survey) ##############################
################################################################################

library(openxlsx)
library(dplyr)
library(ggplot2)
library(ggstar)

total_data = read.xlsx("all_field_data-11-21.xlsx", sheet = "Field_group(all species)11", rowNames = T, colNames = T)
total_data$Sample_ID = rownames(total_data)
total_data$Years = as.factor(total_data$Years)
total_data$Site = factor(total_data$Site, levels = unique(total_data$Site[order(total_data$Latitude)]))

#
total_data_mean = total_data %>% group_by(Site, Years) %>%
  summarise(
    Field_Di = mean(Fungal_field_Di, na.rm = TRUE),
    Field_Di_se = sd(Fungal_field_Di, na.rm = TRUE) / sqrt(n()),
    Green_Di = mean(Fungal_green_Di, na.rm = TRUE),
    Green_Di_se = sd(Fungal_green_Di, na.rm = TRUE) / sqrt(n())) %>%
  as.data.frame()

total_data2 =  merge(total_data, total_data_mean, by = c("Years", "Site"))

ggplot() + 
  geom_star(data = subset(total_data, Years == "2018"), aes(x =Fungal_green_Di , y = Fungal_field_Di, color = Years, fill = Years, starshape = Site), 
            size = 1.2, alpha = 0.6, show.legend = F) + 
  scale_starshape_manual(values = c(13,11,23,15,1,28)) +
  geom_segment(subset(total_data2, Years == "2018"),mapping = aes(x = Green_Di, y = Field_Di, xend = Fungal_green_Di, yend = Fungal_field_Di),
               linetype = 1, linewidth = 0.2, alpha = 0.2, color = "#A38E89")+
  geom_errorbar(data = subset(total_data_mean, Years == "2018"),mapping = aes(x = Green_Di,ymax = Field_Di+Field_Di_se, ymin=Field_Di-Field_Di_se),width=0.008,alpha = 1, color = "black")+
  geom_errorbarh(data = subset(total_data_mean, Years == "2018"),mapping = aes(y = Field_Di,xmax=Green_Di+Green_Di_se,xmin=Green_Di-Green_Di_se),height=0.008,alpha = 1, color = "black")+
  geom_star(data = subset(total_data_mean, Years == "2018"),mapping = aes(x = Green_Di, y = Field_Di, fill = Years, starshape = Site),size=2.5, color = "black")+
  geom_abline(intercept=0,slope=1 , linetype = 2, color = "black")+
  mytheme + #theme(legend.position = "right") + 
  labs(y = "Fungal composition distinctiveness\nestimated in situ observations", x = NULL, tag = "a") + 
  scale_color_manual(values = c("2018" = "#A38E89", "2020" = "#FD7541", "2021" = "#40B0A6")) + 
  scale_fill_manual(values = c("2018" = "#A38E89", "2020" = "#FD7541", "2021" = "#40B0A6")) +
  scale_x_continuous(labels = scales::label_comma(accuracy =0.01), limits = c(0.59,0.94)) +
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01), limits = c(0.59,0.94)) -> p3a;p3a



ggplot() + 
  geom_star(data = subset(total_data, Years == "2020"), aes(x =Fungal_green_Di , y = Fungal_field_Di, color = Years, fill = Years, starshape = Site), 
            size = 1.2, alpha = 0.6, show.legend = F) + 
  scale_starshape_manual(values = c(13,11,23,15,1,28)) +
  geom_segment(subset(total_data2, Years == "2020"),mapping = aes(x = Green_Di, y = Field_Di, xend = Fungal_green_Di, yend = Fungal_field_Di),
               linetype = 1, linewidth = 0.2, alpha = 0.2, color = "#FD7541")+
  geom_errorbar(data = subset(total_data_mean, Years == "2020"),mapping = aes(x = Green_Di,ymax = Field_Di+Field_Di_se, ymin=Field_Di-Field_Di_se),width=0.008,alpha = 1, color = "black")+
  geom_errorbarh(data = subset(total_data_mean, Years == "2020"),mapping = aes(y = Field_Di,xmax=Green_Di+Green_Di_se,xmin=Green_Di-Green_Di_se),height=0.008,alpha = 1, color = "black")+
  geom_star(data = subset(total_data_mean, Years == "2020"),mapping = aes(x = Green_Di, y = Field_Di, fill = Years, starshape = Site),size=2.5, color = "black")+
  geom_abline(intercept=0,slope=1 , linetype = 2, color = "black")+
  mytheme + #theme(legend.position = "right") + 
  labs(y = NULL,
       x = "Fungal composition distinctiveness\nestimated in pot experiment", tag = "b") + 
  scale_color_manual(values = c("2018" = "#A38E89", "2020" = "#FD7541", "2021" = "#40B0A6")) + 
  scale_fill_manual(values = c("2018" = "#A38E89", "2020" = "#FD7541", "2021" = "#40B0A6")) + 
  scale_x_continuous(labels = scales::label_comma(accuracy =0.01), limits = c(0.59,0.94)) +
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01), limits = c(0.59,0.94)) -> p3b;p3b


ggplot() + 
  geom_star(data = subset(total_data, Years == "2021"), aes(x =Fungal_green_Di , y = Fungal_field_Di, color = Years, fill = Years, starshape = Site), 
            size = 1.2, alpha = 0.6, show.legend = F) + 
  scale_starshape_manual(values = c(13,11,23,15,1,28)) +
  geom_segment(subset(total_data2, Years == "2021"),mapping = aes(x = Green_Di, y = Field_Di, xend = Fungal_green_Di, yend = Fungal_field_Di),
               linetype = 1, linewidth = 0.2, alpha = 0.2, color = "#40B0A6")+
  geom_errorbar(data = subset(total_data_mean, Years == "2021"),mapping = aes(x = Green_Di,ymax = Field_Di+Field_Di_se, ymin=Field_Di-Field_Di_se),width=0.008,alpha = 1, color = "black")+
  geom_errorbarh(data = subset(total_data_mean, Years == "2021"),mapping = aes(y = Field_Di,xmax=Green_Di+Green_Di_se,xmin=Green_Di-Green_Di_se),height=0.008,alpha = 1, color = "black")+
  geom_star(data = subset(total_data_mean, Years == "2021"),mapping = aes(x = Green_Di, y = Field_Di, fill = Years, starshape = Site),size=2.5, color = "black")+
  geom_abline(intercept=0,slope=1 , linetype = 2, color = "black")+
  mytheme + #theme(legend.position = "right") + 
  labs(x = NULL, y = NULL, tag = "c")+ 
  scale_color_manual(values = c("2018" = "#A38E89", "2020" = "#FD7541", "2021" = "#40B0A6")) + 
  scale_fill_manual(values = c("2018" = "#A38E89", "2020" = "#FD7541", "2021" = "#40B0A6")) + 
  scale_x_continuous(labels = scales::label_comma(accuracy =0.01), limits = c(0.59,0.94)) +
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01), limits = c(0.59,0.94)) -> p3c;p3c



#### Environmental effects
total_data$Effect_size = log(total_data$Fungal_field_Di/total_data$Fungal_green_Di)

### T.test
Years = unique(total_data$Years)
Latitude = unique(total_data$Latitude)
final_result_t = NULL
for (Y in Years) {
  for (L in Latitude) {
    group_di = subset(total_data, Years == Y & Latitude == L)
    X = group_di$Effect_size
    # 
    result_t = (t.test(X, mu = 0))
    result_t_test = data.frame(Years = Y, Latitude = L, t_value = result_t$statistic,
                               p_value = result_t$p.value)
    final_result_t = rbind(final_result_t, result_t_test)
    rownames(final_result_t) = NULL
  }
}
final_result_t$p_value = round(final_result_t$p_value, 3)

mean_size = total_data %>% group_by(Latitude, Years) %>%
  summarise(
    Effect_mean = mean(Effect_size, na.rm = TRUE),
    Effect_se = sd(Effect_size, na.rm = TRUE) / sqrt(n())) %>%
  left_join(final_result_t, by = c("Latitude", "Years"))
mean_size$Years = as.factor(mean_size$Years)
mean_size$p_value = round(mean_size$p_value, 3)
mean_size$sig <- ifelse(mean_size$p_value>0.05, 0, 1)
mean_size$sig = as.factor(mean_size$sig)
total_data$Years = as.factor(total_data$Years)

summary(mean_size$Effect_mean)

library(BestFitM)
library(ggtrendline)

bestFitM2(data= mean_size, x= "Latitude", y = "Effect_mean")
ggtrendline(mean_size$Latitude, mean_size$Effect_mean, model = "line3P") + 
  geom_point(mean_size, mapping = aes(x = Latitude, y = Effect_mean))

ggplot() + 
  geom_errorbar(data = subset(mean_size, Years == 2018),mapping = aes(x = Latitude-0.15,ymax = Effect_mean+Effect_se, ymin=Effect_mean-Effect_se, color = Years),width=0.2,alpha = 1, color = "black")+
  geom_errorbar(data = subset(mean_size, Years == 2020),mapping = aes(x = Latitude,ymax = Effect_mean+Effect_se, ymin=Effect_mean-Effect_se, color = Years),width=0.2,alpha = 1, color = "black")+
  geom_errorbar(data = subset(mean_size, Years == 2021),mapping = aes(x = Latitude+0.15,ymax = Effect_mean+Effect_se, ymin=Effect_mean-Effect_se, color = Years),width=0.2,alpha = 1, color = "black")+
  geom_star(data = subset(mean_size, Years == 2018 & Latitude != 23.1 & Latitude != 34.6), mapping = aes(x = Latitude-0.15, y = Effect_mean, starshape = factor(Latitude), fill = Years),size = 2.8, show.legend = F) + 
  geom_star(data = subset(mean_size, Years == 2018 & Latitude == 23.1), mapping = aes(x = Latitude-0.15, y = Effect_mean),size = 2.8, show.legend = F, color = "#A38E89", fill = "white", starshape = 13) + 
  geom_star(data = subset(mean_size, Years == 2018 & Latitude == 34.6), mapping = aes(x = Latitude-0.15, y = Effect_mean),size = 2.8, show.legend = F, color = "#A38E89", fill = "white", starshape = 1) + 
  geom_star(data = subset(mean_size, Years == 2020 & Latitude != 30.5 & Latitude != 34.6), mapping = aes(x = Latitude, y = Effect_mean, starshape = factor(Latitude), fill = Years),size = 2.8, show.legend = F) + 
  geom_star(data = subset(mean_size, Years == 2020 & Latitude == 30.5), mapping = aes(x = Latitude, y = Effect_mean),size = 2.8, show.legend = F, color = "#FD7541", fill = "white", starshape = 15) + 
  geom_star(data = subset(mean_size, Years == 2020 & Latitude == 34.6), mapping = aes(x = Latitude, y = Effect_mean),size = 2.8, show.legend = F, color = "#FD7541", fill = "white", starshape = 1) +
  geom_star(data = subset(mean_size, Years == 2021 & Latitude != 36.2), mapping = aes(x = Latitude+0.15, y = Effect_mean, starshape = factor(Latitude), fill = Years),size = 2.8, show.legend = F) + 
  geom_star(data = subset(mean_size, Years == 2021 & Latitude == 36.2), mapping = aes(x = Latitude+0.15, y = Effect_mean),size = 2.8, show.legend = F, color = "#40B0A6", fill = "white", starshape = 28) + 
  geom_smooth(data = total_data, mapping = aes(x = Latitude, y = Effect_size),
              method = "lm", formula = y ~ poly(x, 2), se = F, color = "black")  +
  scale_starshape_manual(values = c("23.1" = 13, "25.2" = 11, "27.9" = 23, "30.5" = 15, "34.6" = 1, "36.2" = 28)) +
  scale_color_manual(values = c("2018" = "#A38E89", "2020" = "#FD7541", "2021" = "#40B0A6")) + 
  scale_fill_manual(values = c("2018" = "#A38E89", "2020" = "#FD7541", "2021" = "#40B0A6")) + 
  geom_hline(yintercept = 0, linetype = 2) +
  mytheme + 
  #annotate("segment", x = 33, xend = 33, y = 0, yend = 0.2, arrow = arrow(angle = 25, length = unit(0.1,  "in")))+ 
  theme(legend.position = c(0.35,0.25)) + 
  labs(x = "Latitude (North degrees)", y = "Effect size", tag = "b")


### Notice that,
### For more picture details, we have further adjusted it in Adobe illustrator.


