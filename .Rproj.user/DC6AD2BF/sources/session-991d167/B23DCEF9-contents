################################################################i
# Title: Site map for Barkley Sound 
#
# Author: Lydia Walton
# Last updated: 13-Dec-2024
################################################################i

# Setup ----
rm(list = ls(all=T)) #Clear environment
#Set working directory
library(here)
#Data manipulation
library(dplyr)
library(tidyr)
#Data visualization
library(ggplot2)
library(kableExtra)
library(patchwork)
#Making the site map
library(maps) 
library(mapdata)
library(tidyverse)
library(PBSmapping)
library(cowplot)
library(mapproj)
library(sf)
library(ggspatial)

#Create theme for plotting

#Lydia's theme
LW_theme <- theme_classic() +
  theme(axis.title.x = element_text(vjust = -1, size = 12),
        axis.title.y = element_text(vjust = 2, size = 12),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 12),
        panel.background = element_rect(fill = NA, color = "black"))


###################i
# SITE MAP ----
###################i

## Van Isle ----

#Use this data for the Van Isle map 
#data(nepacLLhigh)
VanIsle.sf <- read_sf("2023 data/bc-coast.shp")

#Figure out the coordinates 
ggplot() +
  geom_sf(data = VanIsle.sf, fill = "grey75") +
  coord_sf(xlim = c(-129, -122), ylim = c(47.5, 51.5)) +
  LW_theme

# crop map to show the west coast of Canada
VanIsle.Inset <- ggplot() +
  geom_sf(data = VanIsle.sf, fill = "grey75", color = "white") +
  coord_sf(xlim = c(-128.5, -122), ylim = c(48, 51.5)) +
  LW_theme +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.margin=grid::unit(c(0,0,0,0), "mm"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  #Add text for Canada and Vancouver Island
  annotate("text", x = -123, y = 51, label = "Canada", 
           size = 5, fontface = "bold") +
  annotate("text", x = -126.3, y = 50.1, label = "Vancouver \nIsland", 
           size = 4, fontface = "bold") +
  annotate("text", x = -125.25, y = 48.5, label = "Barkley \nSound", 
           size = 3, fontface = "bold") +
  #Add box where Site map will be
  annotate("rect",xmin = -124.9, xmax = -125.55, ymin = 48.75, ymax = 49.1, 
           col = "red", fill = NA, size = 0.8)

VanIsle.Inset


## Barkely sound experiment and field sites ----

#Read in sites
POsites <- read.csv("2023 Data/Summer2023_SiteMap.csv")

#Filter for Sitemap (sites to be included in this map)
POsites.field <- POsites %>% 
  filter(Map == "Sitemap")

#Make site description an ordered factor for plotting
POsites.field$Site.description <- factor(POsites.field$Site.description,
                                   ordered = TRUE,
                                   levels = c("Collections","Experiments", "Moribundity surveys"))

#Hakai basemap
hakai.basemap <- read_sf("2023 Data/COAST_TEST2.shp")
st_crs(hakai.basemap)


#Crop the map for Barkley Sound 
barkley.sound.field <- st_crop(hakai.basemap,
                               c(xmin = -125.08, xmax = -125.18, ymin = 48.81, ymax = 48.853))

#Plot the map to see which area it selects
#ggplot(barkley.sound.temp) +
#  geom_sf(fill = "grey75", color = "white") +
#  annotation_scale(location = "bl", width_hint = 0.4) +
#  coord_sf(expand = FALSE) +
#  theme_void()

###Field sites (MS Figure) ----
SiteMap <- ggplot() + 
  LW_theme +
  geom_sf(data = barkley.sound.field, color= "white", fill='grey75') +
  annotation_scale(location = "bl", width_hint = 0.4) +
  coord_sf(expand = FALSE) +
  labs(x = "Longitude", y = "Latitude") +
#Add site locations
  geom_point(data=POsites.field, 
             aes(x=Lon, y=Lat, 
             color = Site.description, shape = Site.description), 
             size=4) +
  scale_colour_manual(values = c("#0F8B8D", "#041B15", "#9B5094")) +
  scale_shape_manual(values = c(17,18,16)) +
  labs(shape = "Site description", color = "Site description") +
  theme(legend.position = "top",
        #legend.justification = "top",
        axis.title = element_text(size = 12)) +
  annotation_north_arrow(location = "tl", 
                           pad_x = unit(0.4, "in"), pad_y = unit(0.4, "in"),
                           style = north_arrow_nautical(
                           fill = c("grey40", "white"),
                           line_col = "grey20"))
#Add site labels 
SiteMap <- SiteMap +
  geom_label(aes(label= "BMSC", x=-125.1352, y= 48.8368), size=4) +
  geom_label(aes(label= "EB", x=-125.147, y= 48.8352), size=4) +
  geom_label(aes(label= "GN", x=-125.1109, y= 48.8358), size=4) +
  geom_label(aes(label= "SP", x=-125.1296, y= 48.8338), size=4) +
  geom_label(aes(label= "Bamfield Inlet", x=-125.1385, y= 48.8265), size=4) +
  geom_text(aes(label= "Bamfield", x=-125.10, y= 48.845), size=5, fontface = "bold") +
  geom_text(aes(label= "Trevor \nChannel", x=-125.154, y= 48.841), size=5, fontface = "bold")


SiteMap


##### MS Figure ----

#Add van isle inset
#Combine the plots, adding VI map as an inset in the top right corner
Site.map.withinset <- ggdraw() +
  draw_plot(SiteMap) +
  draw_plot(VanIsle.Inset, x = 0.64, y = 0.13, width = 0.35, height = 0.35)

Site.map.withinset

#png("MS Figures/Fig1_Site-map-with-inset.png", width = 11, height = 9, units = "in", res = 600)
#Site.map.withinset
#dev.off()

