library(openxlsx)
seed_data = read.xlsx("Site_information_seed.xlsx", sheet = "site_map", rowNames = F, colNames = T)
head(seed_data)

options(amap_key = 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxx') # Your KEY
library(tidygeocoder)

library(sf)
library(ggplot2)
library(ggspatial)
china_map <- sf::st_read("https://geo.datav.aliyun.com/areas_v3/bound/100000_full.json") 
class(china_map)

###
main_theme = theme(panel.background=element_blank(),
                   panel.grid=element_blank(),
                   axis.line.x=element_line(size=0.5, colour="NA"),
                   axis.line.y=element_line(size=0.5, colour="NA"),
                   axis.ticks=element_line(color="black"),
                   axis.text=element_text(color="black", size=11),
                   legend.position="right",
                   legend.background=element_blank(),
                   legend.key=element_blank(),
                   legend.text= element_text(size=9),
                   legend.title = element_text(size=10),
                   plot.tag = element_text(size = 14, face = "bold"))

ggplot(china_map)+
  geom_sf(data = china_map,fill="grey95",size=1) + 
  xlim(105,122)+ ylim(18,42)+ 
  ggnewscale::new_scale_fill() + 
  #geom_text(data=datasel,aes(x=lon, y=lat-0.6 ,label=province),size=3,colour="black")+
  geom_point(data = seed_site_data, aes(x=Longitude,y=Latitude),
             size=3,alpha=1, shape = 21, color = "black", fill = "#32B0A5") + 
  main_theme+
  annotation_scale(location = "br", style = "ticks",line_width = 0.1,pad_y = unit(0.5, "cm"), text_cex = 0.9) + 
  annotation_north_arrow(location = "tl", which_north = "true", height = unit(1, "cm"),width = unit(1, "cm"),
                         pad_x = unit(0.3, "cm"), pad_y = unit(0.3, "cm"),style = north_arrow_fancy_orienteering) +
  #coord_sf(crs = "+proj=laea +lat_0=40 +lon_0=104")+
  guides(fill = guide_legend(title = "Site",ncol = 1, nrow = 14, override.aes = list(shape = 22, size=5))) +
  theme(text = element_text(size = 13),
        panel.background = element_rect(fill='#FFFFFF', colour='black'),
        axis.line.x = element_line(size=0.5, colour="black"),
        axis.line.x.top =element_line(size=0.5, colour="black"),
        axis.line.y.right = element_line(size=0.5, colour="black"),
        axis.line.y=element_line(size=0.5, colour="black"),
        legend.position = "right",
        panel.grid.major = element_line(color = "white",size=0.2))+
  labs(x='', y='') +
  ggrepel::geom_text_repel(mapping = aes(x=Longitude,y=Latitude,label=`City`), data = seed_site_data,size = 2.8,segment.color = "black", 
                           color = "black",direction = "both",box.padding = 0.6,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 25))
                           

####
library(ggtree)
library(phytools)
tree <-read.newick("iq_tree.NEWICK")
to_drop<-c("Amborella_trichopoda","")
tree <- drop.tip(as.phylo(tree), to_drop) 

ggtree(tree, size = 0.4, color = "black", branch.length = "branch.length", ladderize = F) + 
  geom_tiplab(aes(label = sub("_", " ", label)), size = 2, offset = 0.01,
              fontface = "italic", linetype = "dotted", align = tT) + 
  xlim(0,8)-> p.tree; p.tree

ggtree::rotate(p1, node = 56)

# Add collection point information
seed_site_infor = read.xlsx("Site_information_seed.xlsx", sheet = "sp_site", rowNames = F, colNames = T)
head(seed_site_infor)

library(gt)
library(tidyverse)
library(glue)

gt(seed_site_infor)

### Notice that,
### For more picture details, we have further adjusted it in Adobe illustrator.


