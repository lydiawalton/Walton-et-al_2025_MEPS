#######################################################################################################i
# Title: Investigating the impact of contrasting water and air temperatures on the body temperature of
#        juvenile Pisaster ochraceus when submerged in seawater and emersed in air
#
# Author: Lydia Walton
# Last updated: 26-July-2024
######################################################################################################i

#dev.off() #Clear open plots
rm(list = ls(all=T)) #Clear environment

#Set working directory
library(here)
#Packages for data manipulation
library(dplyr)
library(tidyr)
#Packages for plotting
library(ggplot2)
library(patchwork)
#Packages for statistics and modelling
library(glmmTMB)
library(DHARMa)
library(sjPlot)

# Setup ----

#------------------------------------------Load data files

PO.collection.May2024 <- read.csv("2024 data/PO_Collection_May-2024.csv")
PO.bodytemps.May2024 <- read.csv("2024 data/PO_Body-temperatures_May-2024.csv")

#-----------------------------------------Create theme for plotting

#Lydia's theme
LW_theme <- theme_classic() +
  theme(axis.title.x = element_text(vjust = -1, size = 12),
        axis.title.y = element_text(vjust = 2, size = 11),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        panel.background = element_rect(fill = NA, color = "black"))


## Data manipulation ----

#Subset the data by the variables we want 
PO.collection.May2024.sub <- PO.collection.May2024 %>% #Collection data
  subset(select = c(star.ID, ind.number, water.treatment, air.treatment, disc.diameter,
                    collection.weight))


PO.bodytemps.May2024 <- PO.bodytemps.May2024 %>% 
  subset(select = c(star.ID, ind.number, water.treatment, air.treatment, trial.day,
                    trial.period, trial.block, measurement.type,
                    disc.temp, arm.1, arm.2, arm.3, arm.4, arm.5, ave.star.temp, aquaria.temp,
                    water.temp, air.temp))

#Merge dfs together
bodytemp.master.May2024 <- merge(PO.collection.May2024.sub, PO.bodytemps.May2024,
                                   by = c("star.ID", "ind.number",
                                          "water.treatment", "air.treatment"))

#Make sure all of the variables are the right structure
bodytemp.master.May2024$water.treatment <- as.factor(bodytemp.master.May2024$water.treatment)
bodytemp.master.May2024$air.treatment <- as.factor(bodytemp.master.May2024$air.treatment)

##RESCALE body size for modelling
bodytemp.master.May2024$disc.diameter_z <- arm::rescale(bodytemp.master.May2024$disc.diameter)

bodytemp.master.May2024 <- bodytemp.master.May2024 %>% 
  mutate(Treatment = case_when(
    water.treatment =="15" & air.treatment=="25" ~ "15°C water/25°C air", 
    water.treatment =="15" & air.treatment=="30" ~ "15°C water/30°C air", 
    water.treatment =="20" & air.treatment=="25" ~ "20°C water/25°C air", 
    water.treatment =="20" & air.treatment=="30" ~ "20°C water/30°C air", 
    TRUE ~ "other"
  )) 

#New column for trial "blocks"
bodytemp.master.May2024 <- bodytemp.master.May2024 %>% 
  mutate(trial.block = case_when(
    trial.day < 11 ~ "Adjustment", 
    trial.day > 10 & trial.day < 21 ~ "First half heat stress", 
    trial.day > 20 & trial.day < 31 ~ "Second half heat stress",
    trial.day > 30 ~ "Recovery",
    TRUE ~ "OTHER"
  ))

#Add group column based on odd or even ind.number
bodytemp.master.May2024 <- bodytemp.master.May2024 %>% 
  mutate(group = case_when(
    Treatment == "15°C water/25°C air" & ind.number %% 2 ~ "15C Group 1",
    Treatment == "15°C water/25°C air" & ind.number %% 2 == 0 ~ "15C Group 2",
    Treatment == "15°C water/30°C air" & ind.number %% 2 ~ "15C Group 1",
    Treatment == "15°C water/30°C air" & ind.number %% 2 == 0 ~ "15C Group 2",
    Treatment == "20°C water/25°C air" & ind.number %% 2 ~ "20C Group 1",
    Treatment == "20°C water/25°C air" & ind.number %% 2 == 0 ~ "20C Group 2",
    Treatment == "20°C water/30°C air" & ind.number %% 2 ~ "20C Group 1",
    Treatment == "20°C water/30°C air" & ind.number %% 2 == 0 ~ "20C Group 2", 
    TRUE ~ "other"
  )) %>% 
  drop_na(disc.temp)

# New facet label names for variables
measurement.labs <- c("Emersed", "Submerged")
names(measurement.labs) <- c("emersed", "submerged")

trial.labs <- c("Heat stress", "Recovery")
names(trial.labs) <- c("Heat stress", "Recovery")

## Body temperatures ----

### Statistics ----

# Look at the submerged temperatures first (15C versus 20C)

bodytemps.submerged <- bodytemp.master.May2024 %>% 
  filter(measurement.type == "submerged")

# Correlation test between ave.star.temp and water.temp (actual temperature of the water NOT the set treatment temp)

# Shapiro-Wilk normality test
shapiro.test(bodytemps.submerged$water.temp) #p < 0.001
# Shapiro-Wilk normality test 
shapiro.test(bodytemps.submerged$ave.star.temp) #p < 0.001


bodytemp.submerged.cor <- cor.test(bodytemps.submerged$water.temp, bodytemps.submerged$ave.star.temp,
         method = "kendall", alternative = "greater")

bodytemp.submerged.cor #p < 0.001 Sign. correlation between body temp and water temp

# Look at the emersed temperatures (25C versus 30C)
bodytemps.emersed <- bodytemp.master.May2024 %>% 
  filter(measurement.type == "emersed")

# Correlation test between ave.star.temp and air.temp (actual temperature of the air NOT the set treatment temp)

# Shapiro-Wilk normality test
shapiro.test(bodytemps.emersed$air.temp) #p < 0.001
# Shapiro-Wilk normality test 
shapiro.test(bodytemps.emersed$ave.star.temp) #p < 0.001


bodytemp.emersed.cor <- cor.test(bodytemps.emersed$air.temp, bodytemps.emersed$ave.star.temp,
                                   method = "kendall", alternative = "greater")

bodytemp.emersed.cor #p < 0.001 Sign. correlation between body temp and air temp


### Plotting ----

#Make df for panel labels (A-D)

#New df
bodytemp.plotting <- bodytemp.master.May2024 %>% 
  group_by(water.treatment, air.treatment, Treatment,
           trial.period, trial.block, measurement.type) %>% 
  summarise(mean.BT = mean(ave.star.temp), sd.BT = sd(ave.star.temp),
            water.temp = mean(water.temp),
            air.temp = mean(air.temp))

#Plot submersion temperatures
bodytemp.plotting.submerged <- bodytemp.plotting %>% 
  filter(measurement.type == "submerged")

####PANEL A - heat stress ----
bodytemp.plotting.submerged.heatstress <- bodytemp.plotting.submerged %>% 
  filter(trial.period == "Heat stress")

#plot
bodytemp.plot.PanelA <- bodytemp.plotting.submerged.heatstress %>% 
  ggplot() +
  geom_point(aes(x = Treatment, y = mean.BT, color = Treatment, shape = measurement.type), size = 2) +
  geom_errorbar(aes(ymin = mean.BT-sd.BT, ymax = mean.BT+sd.BT, 
                    x = Treatment, color = Treatment), width=.2) +
  #Add reference lines for actual air and water temps
  geom_abline(intercept = 15, slope=0, linetype = 2, color = "black") +
  geom_abline(intercept = 19.2, slope=0, linetype = 2, color = "black") +
  scale_color_manual(name = "Treatment",
                     values = c("#197BBD", "#05299E" , "#FC814A", "#8D0801")) +
  scale_shape_manual(values = c(12),
                     guide = "none") +
  facet_grid(~trial.period,
             labeller = labeller(trial.period = trial.labs)) +
  scale_y_continuous(limits = c(12, 22), breaks = c(15, 20))+
  LW_theme +
  theme(legend.position = "top",
        axis.title.x = element_blank()) +
  ylab("Body temperature (°C)") +
#Add labels for sea table temperatures
  geom_text(label = "15.0°C\nseawater", x = 4.35, y = 15.7, 
            size = 3.5, color = "black") +
  geom_text(label = "19.2°C\nseawater", x = 4.35, y = 19.9, 
            size = 3.5, color = "black") +
  #Add panel label
  geom_text(label= "A", x = 0.6, y = 22, 
            size = 5, color = "black")

bodytemp.plot.PanelA

####PANEL B - recovery ----
bodytemp.plotting.submerged.recovery <- bodytemp.plotting.submerged %>% 
  filter(trial.period == "Recovery")

#plot
bodytemp.plot.PanelB <- bodytemp.plotting.submerged.recovery%>% 
  ggplot() +
  geom_point(aes(x = Treatment, y = mean.BT, color = Treatment, shape = measurement.type), size = 2) +
  geom_errorbar(aes(ymin = mean.BT-sd.BT, ymax = mean.BT+sd.BT, 
                    x = Treatment, color = Treatment), width=.2) +
  #Add reference lines for actual air and water temps
  geom_abline(intercept = 14.7, slope=0, linetype = 2, color = "black") +
  scale_color_manual(name = "Treatment",
                     values = c("#197BBD", "#05299E" , "#FC814A", "#8D0801")) +
  scale_shape_manual(values = c(12),
                     guide = "none") +
  facet_grid(measurement.type~trial.period,
             labeller = labeller(measurement.type = measurement.labs,
                                 trial.period = trial.labs)) +
  scale_y_continuous(limits = c(12, 22), breaks = c(15, 20))+
  LW_theme +
  theme(legend.position = "top",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  ylab("Body temperature (°C)") +
  #Add labels for sea table temperatures
  geom_text(label = "14.7°C\nseawater", x = 4.35, y = 15.4, 
            size = 3.5, color = "black") +
  ##Add recovery period shading
  annotate("rect", fill = "grey", alpha = 0.3, 
           xmin = -Inf, xmax = Inf,
           ymin = -Inf, ymax = Inf) +
  #Add panel label
  geom_text(label= "B", x = 0.6, y = 22, 
            size = 5, color = "black")

bodytemp.plot.PanelB

#Plot emersion temperatures
bodytemp.plotting.emersed <- bodytemp.plotting %>% 
  filter(measurement.type == "emersed")

#### PANEL C - heat stress ----
bodytemp.plotting.emersed.heatstress <- bodytemp.plotting.emersed %>% 
  filter(trial.period == "Heat stress")

#plot
bodytemp.plot.PanelC <- bodytemp.plotting.emersed.heatstress %>% 
  ggplot() +
  geom_point(aes(x = Treatment, y = mean.BT, color = Treatment, shape = measurement.type), size = 2) +
  geom_errorbar(aes(ymin = mean.BT-sd.BT, ymax = mean.BT+sd.BT, 
                    x = Treatment, color = Treatment), width=.2) +
  #Add reference lines for actual air and water temps
  geom_abline(intercept = 24.4, slope=0, linetype = 2, color = "black") +
  geom_abline(intercept = 28.9, slope=0, linetype = 2, color = "black") +
  scale_color_manual(name = "Treatment",
                    values = c("#197BBD", "#05299E" , "#FC814A", "#8D0801")) +
  scale_shape_manual(values = c(9),
                     guide = "none") +
  facet_grid(~trial.period,
             labeller = labeller(trial.period = trial.labs)) +
  scale_y_continuous(limits = c(17, 30), breaks = c(20, 25, 30))+
  LW_theme +
  theme(legend.position = "top") +
  ylab("Body temperature (°C)") +
  #Add labels for sea table temperatures
  geom_text(label = "24.4°C\nair", x = 4.4, y = 25.1, 
            size = 3.5, color = "black") +
  geom_text(label = "28.9°C\nair", x = 4.4, y = 29.6, 
            size = 3.5, color = "black") +
  #Add panel label
  geom_text(label= "C", x = 0.6, y = 30, 
            size = 5, color = "black")

  
bodytemp.plot.PanelC

#### PANEL D - recovery ----
bodytemp.plotting.emersed.recovery <- bodytemp.plotting.emersed %>% 
  filter(trial.period == "Recovery")

#plot
bodytemp.plot.PanelD <- bodytemp.plotting.emersed.recovery %>% 
  ggplot() +
  geom_point(aes(x = Treatment, y = mean.BT, color = Treatment, shape = measurement.type), size = 2) +
  geom_errorbar(aes(ymin = mean.BT-sd.BT, ymax = mean.BT+sd.BT, 
                    x = Treatment, color = Treatment), width=.2) +
  #Add reference lines for actual air and water temps
  geom_abline(intercept = 20, slope=0, linetype = 2, color = "black") +
  scale_color_manual(name = "Treatment",
                     values = c("#197BBD", "#05299E" , "#FC814A", "#8D0801")) +
  scale_shape_manual(values = c(9),
                     guide = "none") +
  facet_grid(measurement.type~trial.period,
             labeller = labeller(measurement.type = measurement.labs,
                                 trial.period = trial.labs)) +
  scale_y_continuous(limits = c(17, 30), breaks = c(20, 25, 30))+
  LW_theme +
  theme(legend.position = "top",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  #Add labels for sea table temperatures
  geom_text(label = "20.0°C\nair", x = 4.4, y = 20.8, 
            size = 3.5, color = "black") +
  ##Add recovery period shading
  annotate("rect", fill = "grey", alpha = 0.3, 
           xmin = -Inf, xmax = Inf,
           ymin = -Inf, ymax = Inf)  +
  #Add panel label
  geom_text(label= "D", x = 0.6, y = 30, 
            size = 5, color = "black")


bodytemp.plot.PanelD

#Plot for MS ----
MS.bodytemp.plot <- bodytemp.plot.PanelA + bodytemp.plot.PanelB + bodytemp.plot.PanelC + bodytemp.plot.PanelD +
  plot_layout(nrow = 2, ncol = 2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top") 

MS.bodytemp.plot

#png("MS Figures/Fig6.Body-temperature_May2024.png", width = 12, height = 10, units = "in", res = 600)
#MS.bodytemp.plot
#dev.off()


## Table ----

bodytemp.table <- bodytemp.plotting %>% 
  subset(select = c(Treatment, trial.period, measurement.type, mean.BT, sd.BT, water.temp, air.temp))

bodytemp.table









