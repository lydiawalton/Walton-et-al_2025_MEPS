################################################################i
# Title: Air and Sea Surface Temperatures for Barkley Sound
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


##########################################i
# Air temperature data ----
#
## From Eagle Bay, Bamfield 
## (Collected by the Baum Lab at UVic)
##########################################i

##Import data
airtemps.EB <- read.csv("SuppMat data/ScottsBay_TempData_Aug2021-Jul2023_BaumLab.csv")

##Split date and time into separate columns
airtemps.EB <- separate(airtemps.EB, DateTime, c("Date", "Time"), sep = " ")

#Quick plot of the data
EB.tempseries <- airtemps.EB %>% 
  ggplot(aes(x = Date, y = Temperature)) +
  geom_line() +
  geom_point() +
  LW_theme +
  labs(x = "Date (yyyy-mm-dd)", y = "Temperature (°C)")+
  #scale_x_continuous(breaks = seq(0, 16, by = 2)) +
  scale_y_continuous(breaks = seq(0, 35, by = 2)) 

EB.tempseries
dev.off()

### Data manipulation ----

# Split date column into 3 (day, month, year)
airtemps.EB.months <- airtemps.EB %>%
  mutate(Year = lubridate::year(Date), 
         Month = lubridate::month(Date), 
         Day = lubridate::day(Date))

#Change month to a factor
airtemps.EB.months$Month <- as.factor(airtemps.EB.months$Month)

#Change year to a factor
airtemps.EB.months$Year <- as.factor(airtemps.EB.months$Year)

##Create a column with mean, min, and max values of temp 

#-----Per Year
airtemps.means.perYear <- airtemps.EB.months %>% 
  group_by(Year) %>% 
  summarise(meanTemp = mean(Temperature),
            sdTemp = sd(Temperature),
            minTemp = min(Temperature),
            maxTemp = max(Temperature))

#-----Per Month
airtemps.means.perMonth <- airtemps.EB.months %>% 
  group_by(Year, Month) %>% 
  summarise(meanTemp = mean(Temperature),
            sdTemp = sd(Temperature),
            minTemp = min(Temperature),
            maxTemp = max(Temperature))

#-----Per Day
airtemps.means.perDay <- airtemps.EB.months %>% 
  group_by(Year, Month, Day) %>% 
  summarise(meanTemp = mean(Temperature),
            sdTemp = sd(Temperature),
            minTemp = min(Temperature),
            maxTemp = max(Temperature))


### Summary tables of air temperature at Eagle Bay ----

#-------------Per Year
airtable.perYear <- airtemps.means.perYear %>%
  kable(caption = "Table 1. Mean, max and minumum temperatures at Scott's Bay shown per year (2021 and 2022).", 
        #Change the column names
        col.names = c("Year", "Mean temp (°C)", "sd temp (°C)", "min temp (°C)", "max temp (°C)"), 
        align = "c") %>%  #Align the data in the centre of the column
  row_spec(row = 0, bold = TRUE) %>% #Make the column headers bold
  kable_classic(full_width = F, html_font = "Arial") 

airtable.perYear

#-------------Per Month
airtable.perMonth <- airtemps.means.perMonth %>% 
  kable(caption = "Table 1. Mean, max and minumum temperatures at Scott's Bay shown per year (2021 and 2022).", 
        #Change the column names
        col.names = c("Year", "Month", 
                      "Mean temp (°C)", "sd temp (°C)", "min temp (°C)", "max temp (°C)"), 
        align = "c") %>%  #Align the data in the centre of the column
  row_spec(row = 0, bold = TRUE) %>% #Make the column headers bold
  kable_classic(full_width = F, html_font = "Arial") 

airtable.perMonth

#-------------Summer (May-August)
airtable.Summer <- airtemps.means.perMonth %>% 
  filter(Month == "5" | Month == "6" | Month == "7" | Month == "8" ) %>% 
  kable(caption = "Table 1. Mean, max and minumum temperatures at Scott's Bay shown per year (2021 and 2022).", 
        #Change the column names
        col.names = c("Year", "Month", 
                      "Mean temp (°C)", "sd temp (°C)", "min temp (°C)", "max temp (°C)"), 
        align = "c") %>%  #Align the data in the centre of the column
  row_spec(row = 0, bold = TRUE) %>% #Make the column headers bold
  kable_classic(full_width = F, html_font = "Arial") 

airtable.Summer


### Time series plots of air temperature at Eagle Bay ----

## Remove 2021 (doesn't include any of the summer months)
airtemps.means.perMonth <- airtemps.means.perMonth %>% 
  filter(Year != 2021)

#------------------Per Month means
EB.tempseries.perMonth <- airtemps.means.perMonth %>% 
  ggplot(aes(x = Month, y = meanTemp, color = Year, group = Year)) +
  geom_point() +
  geom_line() +
  LW_theme +
  labs(x = "Month", y = "Temperature (°C)")+
  #scale_x_continuous(breaks = seq(0, 16, by = 2)) +
  scale_y_continuous(breaks = seq(0, 17, by = 1)) 

EB.tempseries.perMonth
dev.off()

EB.tempseries.meanmaxmin <- airtemps.means.perMonth %>% 
  ggplot() +
  geom_point(aes(x = Month, y = meanTemp, color = Year, group = Year)) +
  geom_line(aes(x = Month, y = meanTemp, color = Year, group = Year)) +
  geom_ribbon(aes(x=Month, ymin = minTemp, ymax = maxTemp, 
                  color = Year, group = Year), alpha = 0.2) +
  facet_wrap(~Year) +
  LW_theme +
  scale_y_continuous(breaks = seq(0, 40, by = 2))

EB.tempseries.meanmaxmin
dev.off()

###################################################i
# Sea Surface Temperature data ----
#
## From Amphitrite Point, Barkley Sound
## (Collected by the Amphitrite Point Lighthouse)
###################################################i

##Import data
Amph.dailytempdata <- 
  read.csv("SuppMat data/AmphitritePoint_DailySST_and_Salinity_1934-2023.csv")


## Data manipulation ----

# Split date column into 3 (day, month, year)
Amph.dailytempdata.sep <- Amph.dailytempdata %>%
  mutate(Year = lubridate::year(Date), 
         Month = lubridate::month(Date), 
         Day = lubridate::day(Date))

##Filter with the data we want (=>2022)
Amph.chapter2 <- Amph.dailytempdata.sep %>% 
  filter(Year >= 2022) %>% 
  filter(Temp.C < 900)

#-----Per Year
Amph.means.perYear <- Amph.chapter2 %>% 
  group_by(Year) %>% 
  summarise(meanTemp = mean(Temp.C),
            sdTemp = sd(Temp.C),
            minTemp = min(Temp.C),
            maxTemp = max(Temp.C))

#-----Per Month
Amph.means.perMonth <- Amph.chapter2 %>% 
  group_by(Year, Month) %>% 
  summarise(meanTemp = mean(Temp.C),
            sdTemp = sd(Temp.C),
            minTemp = min(Temp.C),
            maxTemp = max(Temp.C))

#-----Per Day
Amph.means.perDay <- Amph.chapter2 %>% 
  group_by(Year, Month, Day) %>% 
  summarise(meanTemp = mean(Temp.C),
            sdTemp = sd(Temp.C),
            minTemp = min(Temp.C),
            maxTemp = max(Temp.C))

## Summary tables of SST at Amphitrite Point ----

#Change month to a factor
Amph.means.perDay$Month <- as.factor(Amph.means.perDay$Month)
Amph.means.perMonth$Month <- as.factor(Amph.means.perMonth$Month)

#Change year to a factor
Amph.means.perDay$Year <- as.factor(Amph.means.perDay$Year)
Amph.means.perMonth$Year <- as.factor(Amph.means.perMonth$Year)

#-------------Per Year
SSTtable.perYear <- Amph.means.perYear %>%
  kable(caption = "Table 1. Mean, max and minumum temperatures at Amphitrite Point shown per year (2021 and 2022).", 
        #Change the column names
        col.names = c("Year", "Mean temp (°C)", "sd temp (°C)", "min temp (°C)", "max temp (°C)"), 
        align = "c") %>%  #Align the data in the centre of the column
  row_spec(row = 0, bold = TRUE) %>% #Make the column headers bold
  kable_classic(full_width = F, html_font = "Arial") 

SSTtable.perYear

#-------------Per Month
SSTtable.perMonth <- Amph.means.perMonth %>% 
  kable(caption = "Table 1. Mean, max and minumum temperatures at Amphitrite Point shown per year (2021 and 2022).", 
        #Change the column names
        col.names = c("Year", "Month", 
                      "Mean temp (°C)", "sd temp (°C)", "min temp (°C)", "max temp (°C)"), 
        align = "c") %>%  #Align the data in the centre of the column
  row_spec(row = 0, bold = TRUE) %>% #Make the column headers bold
  kable_classic(full_width = F, html_font = "Arial") 

SSTtable.perMonth

#-------------Summer (May-August)
SSTtable.Summer <- Amph.means.perMonth %>% 
  filter(Month == "5" | Month == "6" | Month == "7" | Month == "8" ) %>% 
  kable(caption = "Table 1. Mean, max and minumum temperatures at Amphitrite Point shown per year (2021 and 2022).", 
        #Change the column names
        col.names = c("Year", "Month", 
                      "Mean temp (°C)", "sd temp (°C)", "min temp (°C)", "max temp (°C)"), 
        align = "c") %>%  #Align the data in the centre of the column
  row_spec(row = 0, bold = TRUE) %>% #Make the column headers bold
  kable_classic(full_width = F, html_font = "Arial") 

SSTtable.Summer


## Time series plots of SST at Amphitrite Point ----

#------------------Per Month means
AP.tempseries.perMonth <- Amph.means.perMonth %>% 
  ggplot(aes(x = Month, y = meanTemp, color = Year, group = Year)) +
  geom_point() +
  geom_line() +
  LW_theme +
  labs(x = "Month", y = "Temperature (°C)")+
  #scale_x_continuous(breaks = seq(0, 16, by = 2)) +
  scale_y_continuous(breaks = seq(0, 17, by = 1)) 

AP.tempseries.perMonth
dev.off()

AP.tempseries.meanmaxmin <- Amph.means.perMonth %>% 
  ggplot() +
  geom_point(aes(x = Month, y = meanTemp, color = Year, group = Year)) +
  geom_line(aes(x = Month, y = meanTemp, color = Year, group = Year)) +
  geom_ribbon(aes(x=Month, ymin = minTemp, ymax = maxTemp, 
                  color = Year, group = Year), alpha = 0.2) +
  facet_wrap(~Year) +
  LW_theme +
  scale_y_continuous(breaks = seq(0, 40, by = 2))

AP.tempseries.meanmaxmin
dev.off()

#############################################i
# SuppMat Temperature Figure ----
#############################################i

## Panel A ----

#Add column to distinguish air and water temps (Month df)
airtemps.means.perMonth <- airtemps.means.perMonth %>% 
  mutate(measurement.type = "Air")

Amph.means.perMonth <- Amph.means.perMonth %>% 
  mutate(measurement.type = "SST")

#Combine the air and SST data into one df
air.SST.permonth <- rbind(airtemps.means.perMonth, Amph.means.perMonth)

#Add column to distinguish air and water temps (Month df)
airtemps.means.perDay <- airtemps.means.perDay %>% 
  mutate(measurement.type = "Air")

Amph.means.perDay <- Amph.means.perDay %>% 
  mutate(measurement.type = "SST")

#Combine the air and SST data into one df
air.SST.perday <- rbind(airtemps.means.perDay, Amph.means.perDay)

#Filter for just summer months
air.SST.perday.SUMMER <- air.SST.perday %>% 
  filter(Month %in% (5:8)) %>% 
  filter(Year == 2022)

#New labels for facet_grid
month.labels <- c("May", "June", "July", "August")
names(month.labels) <- c("5", "6", "7", "8")

#Temp series across months
tempseries.bymonth <- air.SST.permonth %>% 
  ggplot() +
  geom_point(aes(x = Month, y = meanTemp, color = measurement.type, group = measurement.type),
             size = 1) +
  geom_errorbar(aes(x= Month, ymin = minTemp, ymax = maxTemp, 
                    color = measurement.type, group = measurement.type), 
                linewidth = 0.5) +
  scale_color_manual(values = c("#FFA3AF", "#96C5F7")) +
  facet_wrap(~Year) +
  LW_theme +
  scale_y_continuous(breaks = seq(0, 40, by = 2)) +
  labs(y = "Temperature (°C)", x = "Day", color = "Measurement type") +
  theme(legend.position = "top")

tempseries.bymonth

#Temp series across days (note that SST only has one value per day while air has multiple)
#------------2022 only (parameterizing temperature treatments)
tempseries.byday.maxtemp <- air.SST.perday.SUMMER %>% 
  ggplot() +
  geom_point(aes(x = Day, y = maxTemp, color = measurement.type, group = measurement.type),
             size = 1) +
  geom_line(aes(x = Day, y = maxTemp, color = measurement.type, group = measurement.type)) +
  scale_color_manual(values = c("#FFA3AF", "#96C5F7")) +
  facet_wrap(~Month,
             labeller = labeller(Month = month.labels),
             nrow = 2) +
  LW_theme +
  scale_y_continuous(breaks = seq(4, 32, by = 2), limits = c(4,32)) +
  scale_x_continuous(breaks = seq(0, 31, by = 5)) +
  labs(y = "Temperature (°C)", x = "Day of the month", color = "Measurement type") +
  theme(legend.position = "top")

tempseries.byday.maxtemp


## Panel B ----

#Import site data
POsites <- read.csv("2023 Data/Summer2023_SiteMap.csv")

##New df with temp sites
POsites.temp <- POsites %>% 
  filter(Map == "Tempmap")

### Van Isle Inset ----

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
  annotate("text", x = -125.25, y = 48.45, label = "Barkley \nSound", 
           size = 3, fontface = "bold") +
  #Add box where Site map will be
  annotate("rect",xmin = -124.9, xmax = -125.55, ymin = 48.75, ymax = 49.1, 
           col = "red", fill = NA, size = 0.8)

VanIsle.Inset

### Barkley Sound Map ----

#Hakai basemap
hakai.basemap <- read_sf("2023 Data/COAST_TEST2.shp")
st_crs(hakai.basemap)


#Crop the map for Barkley Sound (This is good distance for temp map)
barkley.sound.temp <- st_crop(hakai.basemap, 
                              c(xmin = -125, xmax = -125.65, ymin = 48.75, ymax = 49.1))


#Make the plot
SiteMap.temperature <- ggplot() + 
  LW_theme +
  geom_sf(data = barkley.sound.temp, color= "white", fill='grey75') +
  annotation_scale(location = "bl", width_hint = 0.4) +
  coord_sf(expand = FALSE) +
  labs(x = "Longitude", y = "Latitude") +
  #Add site locations
  geom_point(data=POsites.temp, 
             aes(x=Lon, y=Lat, 
                 color = Site.description), 
             size=4) +
  scale_colour_manual(values = c("#FFA3AF", "#96C5F7")) +
  labs(shape = "Site description", color = "Site description") +
  theme(legend.position = "right",
        axis.title = element_text(size = 12)) +
  annotation_north_arrow(location = "tl", 
                         pad_x = unit(0.1, "in"), pad_y = unit(0.1, "in"),
                         style = north_arrow_nautical(
                           fill = c("grey40", "white"),
                           line_col = "grey20"))
#Add site labels 
SiteMap.temperature <- SiteMap.temperature +
  geom_label(aes(label= "Amphitrite Point", x=-125.54, y= 48.905), size=5) +
  geom_label(aes(label= "Eagle Bay", x=-125.147, y= 48.817), size=5)

SiteMap.temperature


# SUPP FIGURE ----

#Add van isle inset
#Combine the plots, adding VI map as an inset in the top right corner
Supp.map.withinset <- ggdraw() +
  draw_plot(SiteMap.temperature) +
  draw_plot(VanIsle.Inset, x = 0.54, y = 0.635, width = 0.38, height = 0.38)

Supp.map.withinset

##Site map of air temp and SST locations
Supp.temp.map <- tempseries.byday.maxtemp / Supp.map.withinset +
  plot_layout(widths = 1) +
  plot_annotation(tag_levels = "A") & 
  theme(plot.margin = margin(0,0,1,0, "cm"))

Supp.temp.map

#png("SuppMat Figures/FigS1_Supp-temp-map_2022.png", width = 9, height = 13, units = "in", res = 600)
#Supp.temp.map
#dev.off()



