##############################################################################################i
# Title: Investigating the impact of contrasting water and air temperatures on the feeding
#        rate of juvenile Pisaster ochraceus
#
# Author: Lydia Walton
# Last updated: 13-Dec-2024
############################################################################################i

# Setup ----
#dev.off() #Clear open plots
rm(list = ls(all=T)) #Clear environment

#Set working directory
library(here)
#Packages for data manipulation
library(dplyr)
library(tidyr)
#Packages for statistics and modeling
library(AICcmodavg)
library(glmmTMB)
library(EnvStats)
library(broom)
library(DHARMa)
library(emmeans)
#Packages for data visualization
library(ggplot2)
library(ggpmisc)
library(tibble)
library(quantreg)
library(kableExtra)
library(sjPlot)
library(RColorBrewer)
library(extrafont)
library(ggpubr)
library(patchwork)

#------------------------------------------Load data files

PO.collection.June2023 <- read.csv("2023 Data/Pisaster_Collections_June2023.csv")
PO.heatstress.June2023 <- read.csv("2023 Data/Pisaster_HeatstressPeriod_June2023.csv")
PO.recovery.June2023 <- read.csv("2023 Data/Pisaster_RecoveryPeriod_June2023.csv")

PO.collection.May2024 <- read.csv("2024 data/PO_Collection_May-2024.csv")
PO.adjustment.May2024 <- read.csv("2024 data/PO_Adjustment_May-2024.csv")
PO.heat.stress.May2024 <- read.csv("2024 data/PO_Heatstress_May-2024.csv")
PO.recovery.May2024 <- read.csv("2024 data/PO_Recovery_May-2024.csv")

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


#########################i
# Experiment 1 ----
# (June 2023)
#########################i

## Data cleaning and manipulation ----

#Subset the data by the variables we want 
PO.collection.June2023.sub <- PO.collection.June2023 %>% #Collection data
  subset(select = c(star.ID, arm.length, disc.diameter, collection.weight))


PO.heatstress.June2023.sub <- PO.heatstress.June2023 %>% #Adj day (day9) + heat stress
  subset(select=c(star.ID, water.treatment, air.treatment, trial.period, heat.trial.day, total.mussels.consumed, cum.mussels, trial.type, group))


PO.recovery.June2023.sub <- PO.recovery.June2023 %>% #Recovery
  subset(select=c(star.ID, water.treatment, air.treatment, trial.period, rec.trial.day, rec.mussels.consumed, cum.mussels, trial.type, group))

#Change the column names to match
colnames(PO.heatstress.June2023.sub) <- c("star.ID", "water.treatment", "air.treatment", "trial.period", "day.of.trial", "mussels.consumed","cum.mussels", "trial.type", "group")

colnames(PO.recovery.June2023.sub) <- c("star.ID", "water.treatment", "air.treatment", "trial.period", "day.of.trial", "mussels.consumed","cum.mussels", "trial.type", "group")

#Merge dfs together
feeding.master.June2023 <- rbind(PO.heatstress.June2023.sub, PO.recovery.June2023.sub)

#Add the collection data
feeding.master.June2023 <- merge(feeding.master.June2023, PO.collection.June2023.sub, by = "star.ID")

#Make sure all of the variables are the right structure
feeding.master.June2023$water.treatment <- as.factor(feeding.master.June2023$water.treatment)
feeding.master.June2023$air.treatment <- as.factor(feeding.master.June2023$air.treatment)

#Remove #67 (4 armed star) and the sea stars that died later
feeding.master.June2023 <- feeding.master.June2023 %>% 
  filter(star.ID != "SB_2023_ind67", star.ID != "SB_2023_ind9", star.ID != "SB_2023_ind11", star.ID != "SB_2023_ind18", star.ID != "SB_2023_ind19", star.ID != "SB_2023_ind21", star.ID != "SB_2023_ind23", star.ID != "SB_2023_ind53", star.ID != "SB_2023_ind58")

#Remove NAs
feeding.master.June2023 <- na.omit(feeding.master.June2023)

##RESCALE body size for modelling
feeding.master.June2023$disc.diameter_z <- arm::rescale(feeding.master.June2023$disc.diameter)

feeding.master.June2023 <- feeding.master.June2023 %>% 
  mutate(Treatment = case_when(
    water.treatment =="15" & air.treatment=="20" ~ "15°C water/20°C air", 
    water.treatment =="15" & air.treatment=="25" ~ "15°C water/25°C air", 
    water.treatment =="15" & air.treatment=="30" ~ "15°C water/30°C air", 
    water.treatment =="20" & air.treatment=="20" ~ "20°C water/20°C air", 
    water.treatment =="20" & air.treatment=="25" ~ "20°C water/25°C air", 
    water.treatment =="20" & air.treatment=="30" ~ "20°C water/30°C air", 
    TRUE ~ "other"
  ))


feeding.master.June2023.modelling <- feeding.master.June2023 %>% 
  filter(trial.period != "Adjustment") %>% 
  group_by(star.ID, day.of.trial, trial.type, disc.diameter_z, 
           Treatment, water.treatment, air.treatment, group) %>% 
  summarise(total.number.mussels = sum(mussels.consumed))


## MODELLING ----

### Overall glmm ----
daily.heatrecovery.glmm <- glmmTMB(total.number.mussels ~ trial.type*Treatment
                                   + disc.diameter_z 
                                   + (1|day.of.trial)
                                   + (1|group),
                                   data = feeding.master.June2023.modelling, 
                                   family = poisson)

summary(daily.heatrecovery.glmm) 

daily.heatrecovery.dharma <- simulateResiduals(fittedModel = daily.heatrecovery.glmm, plot = F)

plot(daily.heatrecovery.dharma)

testDispersion(daily.heatrecovery.glmm, alternative = c("two.sided", "greater", "less"),
               plot = T, type = "DHARMa")

tab_model(daily.heatrecovery.glmm,
          show.est = TRUE) #visualizing model outputs


### Heat stress glmm -----

heatstress.modelling <- feeding.master.June2023 %>% 
  filter(trial.type == "Heat stress") %>% 
  group_by(star.ID, disc.diameter_z, water.treatment, air.treatment, group) %>% 
  summarise(total.number.mussels = sum(mussels.consumed))

#Making the GLMM
heatstress.glmm <- glmmTMB(total.number.mussels ~ water.treatment*air.treatment +
                             + disc.diameter_z
                           + (1|group),
                           data = heatstress.modelling, 
                           family = poisson)

summary(heatstress.glmm) 

heatstress.dharma <- simulateResiduals(fittedModel = heatstress.glmm, plot = F)

plot(heatstress.dharma)

tab_model(heatstress.glmm,
          show.est = TRUE)

### 8-day Blocks ----

###############i
# DAY 1 to 8
###############i
heatstress.day1to8.modelling <- feeding.master.June2023 %>% 
  filter(trial.type == "Heat stress") %>% 
  filter(day.of.trial == "1" | day.of.trial == "2" | day.of.trial == "3" | day.of.trial == "4" |
           day.of.trial == "5" | day.of.trial == "6" | day.of.trial == "7" | day.of.trial == "8") %>% 
  group_by(star.ID, disc.diameter_z, water.treatment, air.treatment, group) %>% 
  summarise(total.number.mussels = sum(mussels.consumed))

#Making the GLMM
heatstress.day1to8.glmm <- glmmTMB(total.number.mussels ~ water.treatment*air.treatment +
                                     + disc.diameter_z
                                   + (1|group),
                                   data = heatstress.day1to8.modelling, 
                                   family = poisson)

summary(heatstress.day1to8.glmm) 

heatstress.day1to8.dharma <- simulateResiduals(fittedModel = heatstress.day1to8.glmm, plot = F)

plot(heatstress.day1to8.dharma)

tab_model(heatstress.day1to8.glmm,
          show.est = TRUE)

#Create an object with the estimates and CI 
heatstress.day1to8.estimates <- plot_model(heatstress.day1to8.glmm, type = "int")


##############i
# Day 9 -16
##############i

heatstress.day9to16.modelling <- feeding.master.June2023 %>% 
  filter(trial.type == "Heat stress") %>% 
  filter(day.of.trial == "9" | day.of.trial == "10" | day.of.trial == "11" | day.of.trial == "12" |
           day.of.trial == "13" | day.of.trial == "14" | day.of.trial == "15" | day.of.trial == "16") %>% 
  group_by(star.ID, disc.diameter_z, water.treatment, air.treatment, group) %>% 
  summarise(total.number.mussels = sum(mussels.consumed))

#Making the GLMM
heatstress.day9to16.glmm <- glmmTMB(total.number.mussels ~ water.treatment*air.treatment +
                                      + disc.diameter_z
                                    + (1|group),
                                    data = heatstress.day9to16.modelling, 
                                    family = poisson)

summary(heatstress.day9to16.glmm) 

heatstress.day9to16.dharma <- simulateResiduals(fittedModel = heatstress.day9to16.glmm, plot = F)

plot(heatstress.day9to16.dharma)

testDispersion(heatstress.day9to16.dharma, alternative = c("two.sided", "greater", "less"),
               plot = T, type = "DHARMa")

tab_model(heatstress.day9to16.glmm,
          show.est = TRUE)

#Create an object with the estimates and CI 
heatstress.day9to16.estimates <- plot_model(heatstress.day9to16.glmm, type = "int")


### Recovery ----

##############i
# Day 17-24
##############i

recovery.modelling <- feeding.master.June2023 %>% 
  filter(trial.type == "Recovery") %>% 
  group_by(star.ID, disc.diameter_z, water.treatment, air.treatment, group) %>% 
  summarise(total.number.mussels = sum(mussels.consumed))

#Making the GLMM
recovery.glmm <- glmmTMB(total.number.mussels ~ water.treatment*air.treatment +
                           + disc.diameter_z
                         + (1|group),
                         data = recovery.modelling, 
                         family = poisson)

summary(recovery.glmm) 

recovery.dharma <- simulateResiduals(fittedModel = recovery.glmm, plot = F)

plot(recovery.dharma)

testDispersion(recovery.glmm, alternative = c("two.sided", "greater", "less"),
               plot = T, type = "DHARMa")

tab_model(recovery.glmm,
          show.est = TRUE)

#Create an object with the estimates and CI 
recovery.estimates <- plot_model(recovery.glmm, type = "int")

## PLOTTING ----

##Set up the data for plotting
heatrecovery_plotting.2023 <- feeding.master.June2023 %>% 
  filter(trial.type != "Adjustment") %>% 
  mutate(trial.block = case_when(
    day.of.trial %in% 1:8 ~ "Day 1-8",
    day.of.trial %in% 9:16 ~ "Day 9-16",
    day.of.trial %in% 17:24 ~ "Day 17-24",
    TRUE ~ "other"
  )) %>% 
  mutate(trial.block.byname = case_when(
    day.of.trial %in% 1:8 ~ "First half heat stress",
    day.of.trial %in% 9:16 ~ "Second half heat stress",
    day.of.trial %in% 17:24 ~ "Recovery",
    TRUE ~ "other"
  )) %>% 
  group_by(star.ID, trial.block, trial.block.byname, trial.type, water.treatment, air.treatment) %>% 
  summarise(total.number.mussels = sum(mussels.consumed)) 


##Make trial block an ordered factor
heatrecovery_plotting.2023$trial.block <- factor(heatrecovery_plotting.2023$trial.block,
                                            ordered = TRUE,
                                            levels = c("Day 1-8", "Day 9-16", "Day 17-24"))
##Make trial block an ordered factor
heatrecovery_plotting.2023$trial.block.byname <- factor(heatrecovery_plotting.2023$trial.block.byname,
                                                   ordered = TRUE,
                                                   levels = c("First half heat stress", 
                                                              "Second half heat stress", 
                                                              "Recovery"))


##Add a treatment column
heatrecovery_plotting.2023 <- heatrecovery_plotting.2023 %>% 
  mutate(Treatment = case_when(
    water.treatment =="15" & air.treatment=="20" ~ "15°C water/20°C air", 
    water.treatment =="15" & air.treatment=="25" ~ "15°C water/25°C air", 
    water.treatment =="15" & air.treatment=="30" ~ "15°C water/30°C air", 
    water.treatment =="20" & air.treatment=="20" ~ "20°C water/20°C air", 
    water.treatment =="20" & air.treatment=="25" ~ "20°C water/25°C air", 
    water.treatment =="20" & air.treatment=="30" ~ "20°C water/30°C air", 
    TRUE ~ "other"
  )) 


###############################i
### Panel A ----
### Total mussel consumption
###############################i

#Load the model estimate data
PO_modelestimates_byBlock.2023 <- read.csv("2023 Data/Pisaster_heat-recovery_model-estimates_by-Trial.Block.csv")

#Make appropriate columns factors
PO_modelestimates_byBlock.2023$trial.block <- factor(PO_modelestimates_byBlock.2023$trial.block,
                                                ordered = TRUE,
                                                levels = c("Day 1-8", "Day 9-16", "Day 17-24"))

#Add new column with trial block by name
PO_modelestimates_byBlock.2023 <- PO_modelestimates_byBlock.2023 %>% 
  mutate(trial.block.byname = case_when(
    trial.block == "Day 1-8" ~ "First half heat stress", 
    trial.block == "Day 9-16" ~ "Second half heat stress",
    trial.block == "Day 17-24" ~ "Recovery",
    TRUE ~ "other"
  ))

##Make trial block an ordered factor
PO_modelestimates_byBlock.2023$trial.block.byname <- factor(PO_modelestimates_byBlock.2023$trial.block.byname,
                                                       ordered = TRUE,
                                                       levels = c("First half heat stress", 
                                                                  "Second half heat stress", 
                                                                  "Recovery"))

PO_modelestimates_byBlock.2023$water.treatment <- as.factor(PO_modelestimates_byBlock.2023$water.treatment)
PO_modelestimates_byBlock.2023$air.treatment <- as.factor(PO_modelestimates_byBlock.2023$air.treatment)

#Add treatment column
PO_modelestimates_byBlock.2023 <- PO_modelestimates_byBlock.2023 %>% 
  mutate(Treatment = case_when(
    water.treatment =="15" & air.treatment=="20" ~ "15°C water/20°C air", 
    water.treatment =="15" & air.treatment=="25" ~ "15°C water/25°C air", 
    water.treatment =="15" & air.treatment=="30" ~ "15°C water/30°C air", 
    water.treatment =="20" & air.treatment=="20" ~ "20°C water/20°C air", 
    water.treatment =="20" & air.treatment=="25" ~ "20°C water/25°C air", 
    water.treatment =="20" & air.treatment=="30" ~ "20°C water/30°C air", 
    TRUE ~ "other"
  ))


#Make the baseplot with the boxes
Overall_PLOT_PanelA <- heatrecovery_plotting.2023 %>% 
  ggplot() +
  geom_point(aes(x = Treatment, y = total.number.mussels, 
                 shape = trial.block.byname, color = Treatment),
             position = position_jitterdodge(dodge.width = 1),
             alpha = 0.5) +
  scale_color_manual(name = "Treatment (°C)", 
                     values = c("#9DC7C8", "#197BBD", "#05299E" , "#FC6DAB", "#FC814A", "#8D0801")) +
  scale_shape_manual(name = "Trial block",
                     values = c(17, 25, 8)) +
  LW_theme +
  theme(legend.position = "top",
        legend.box = "vertical",
        legend.spacing.y = unit(-0.1, "cm")) +
  labs(x = "", y = "Total # of mussels consumed per individual") +
  scale_y_continuous(breaks = c(0,1,2,3,4,5,6,7,8), limits = c(-0.6,8))

##ADDING THE ESTIMATES

Overall_PLOT_PanelA <- Overall_PLOT_PanelA +
  geom_point(data = PO_modelestimates_byBlock.2023,
             aes(x = Treatment, y = estimate, group = trial.block.byname, 
                 color = Treatment),
             position = position_dodge(width = 1), size = 4) +
  scale_color_manual(name = "Treatment (°C)", 
                     values = c("#9DC7C8", "#197BBD", "#05299E" , "#FC6DAB", "#FC814A", "#8D0801")) +
  #Add confidence interval for the model coefficient
  geom_errorbar(data = PO_modelestimates_byBlock.2023,
                aes(x = Treatment, y = estimate, group = trial.block.byname,
                    ymin = conf.low, ymax = conf.high, color = Treatment),
                position = position_dodge(width = 1), width = 0.5, linewidth = 1) +
  #Add sample sizes
  geom_text(label = "n=11", x = 1, y = -0.7, size = 4) +
  geom_text(label = "n=12", x = 2, y = -0.7, size = 4) +
  geom_text(label = "n=7", x = 3, y = -0.7, size = 4) +
  geom_text(label = "n=12", x = 4, y = -0.7, size = 4) +
  geom_text(label = "n=12", x = 5, y = -0.7, size = 4) +
  geom_text(label = "n=9", x = 6, y = -0.7, size = 4) +
  
  ##Add recovery period shading
  annotate("rect", fill = "grey", alpha = 0.3, 
           xmin = 1.2, xmax = 1.5,
           ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "grey", alpha = 0.3, 
           xmin = 2.2, xmax = 2.5,
           ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "grey", alpha = 0.3, 
           xmin = 3.2, xmax = 3.5,
           ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "grey", alpha = 0.3, 
           xmin = 4.2, xmax = 4.5,
           ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "grey", alpha = 0.3, 
           xmin = 5.2, xmax = 5.5,
           ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "grey", alpha = 0.3, 
           xmin = 6.2, xmax = 6.5,
           ymin = -Inf, ymax = Inf)

Overall_PLOT_PanelA


#######################i
# Experiment 2 ----
# May 2024
#######################i

## Data cleaning and manipulation ----

#Subset the data by the variables we want 
PO.collection.May2024.sub <- PO.collection.May2024 %>% #Collection data
  subset(select = c(star.ID, ind.number, water.treatment, air.treatment, disc.diameter,
                    collection.weight))


PO.adjustment.May2024.sub <- PO.adjustment.May2024 %>% #Adjustment period (day 1 - 10)
  subset(select=c(star.ID, water.treatment, air.treatment, date, trial.period, trial.day,
                  ST.end.day,total.mussels.consumed, cumulative.mussels, 
                  MR.measured, lactic.acid.sampled))

PO.heat.stress.May2024.sub <- PO.heat.stress.May2024 %>% #Heat stress period (Day 11 - 29)
  subset(select=c(star.ID, water.treatment, air.treatment, date, trial.period, trial.day,
                  ST.end.day,total.mussels.consumed, cumulative.mussels, 
                  MR.measured, lactic.acid.sampled))

PO.recovery.May2024.sub <- PO.recovery.May2024 %>% #Heat stress period (Day 11 - 29)
  subset(select=c(star.ID, water.treatment, air.treatment, date, trial.period, trial.day,
                  ST.end.day,total.mussels.consumed, cumulative.mussels, 
                  MR.measured, lactic.acid.sampled))

#Merge dfs together
feeding.master.May2024 <-rbind(PO.adjustment.May2024.sub, PO.heat.stress.May2024.sub)

feeding.master.May2024 <-rbind(feeding.master.May2024, PO.recovery.May2024.sub)

feeding.master.May2024 <- merge(PO.collection.May2024.sub, feeding.master.May2024,
                                by = c("star.ID", "water.treatment", "air.treatment"))

#Make sure all of the variables are the right structure
feeding.master.May2024$water.treatment <- as.factor(feeding.master.May2024$water.treatment)
feeding.master.May2024$air.treatment <- as.factor(feeding.master.May2024$air.treatment)

##RESCALE body size for modelling
feeding.master.May2024$disc.diameter_z <- arm::rescale(feeding.master.May2024$disc.diameter)

feeding.master.May2024 <- feeding.master.May2024 %>% 
  mutate(Treatment = case_when(
    water.treatment =="15" & air.treatment=="25" ~ "15°C water/25°C air", 
    water.treatment =="15" & air.treatment=="30" ~ "15°C water/30°C air", 
    water.treatment =="20" & air.treatment=="25" ~ "20°C water/25°C air", 
    water.treatment =="20" & air.treatment=="30" ~ "20°C water/30°C air", 
    TRUE ~ "other"
  ))

#New column for trial "blocks"
feeding.master.May2024 <- feeding.master.May2024 %>% 
  mutate(trial.block = case_when(
    trial.day < 11 ~ "Adjustment", 
    trial.day > 10 & trial.day < 21 ~ "First half heat stress", 
    trial.day > 20 & trial.day < 31 ~ "Second half heat stress",
    trial.day > 30 ~ "Recovery",
    TRUE ~ "OTHER"
  ))

#Add group column based on odd or even ind.number
feeding.master.May2024 <- feeding.master.May2024 %>% 
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
  ))


feeding.master.May2024.modelling <- feeding.master.May2024 %>%
  group_by(star.ID, trial.day, trial.period, trial.block, group, disc.diameter_z, 
           Treatment, water.treatment, air.treatment, MR.measured, lactic.acid.sampled) %>% 
  summarise(total.number.mussels = sum(total.mussels.consumed))

#Make a new df to create a column for handling stress (from MR measured)
handling.df <- feeding.master.May2024 %>% 
  filter(trial.day == "9") %>% 
  mutate(MO2.group = case_when(
    MR.measured == "Y" ~ "Metabolic group",
    MR.measured == "N" ~ "Control group",
    TRUE ~ "other"
  )) %>% 
  subset(select = c(star.ID, Treatment, MO2.group))
#Add to modelling df
feeding.master.May2024.modelling <- merge(feeding.master.May2024.modelling, handling.df,
                                          by = c("star.ID", "Treatment"))

#Remove NAs for when animals were fasted (48 h prior and 24 h after to oxygen consumption measurements)
feeding.master.May2024.modelling %>% drop_na(total.number.mussels)

#Remove random NAs from importing csv
feeding.master.May2024.modelling %>% drop_na(star.ID)


## MODELLING ----

##Overall GLMM
overall.glmm <- glmmTMB(total.number.mussels ~ trial.period*Treatment
                        + disc.diameter_z 
                        + (1|trial.day)
                        + (1|group/star.ID)
                        + (1|MO2.group),
                        data = feeding.master.May2024.modelling, 
                        family = poisson)

summary(overall.glmm) 

#With MR measured not as random effect
overall.glmm.MRmeasured <- glmmTMB(total.number.mussels ~ trial.period*Treatment
                                   + disc.diameter_z + MO2.group
                                   + (1|trial.day)
                                   + (1|group/star.ID),
                                   data = feeding.master.May2024.modelling, 
                                   family = poisson)

summary(overall.glmm.MRmeasured) #DID not converge

#With individual nested in MO2 group
#With MR measured not as random effect
overall.glmm.nested <- glmmTMB(total.number.mussels ~ trial.period*Treatment
                               + disc.diameter_z
                               + (1|trial.day)
                               + (1|MR.measured/star.ID),
                               data = feeding.master.May2024.modelling, 
                               family = poisson)

summary(overall.glmm.nested)

#Compare models
AIC(overall.glmm.MRmeasured,overall.glmm, overall.glmm.nested)
##Including MO2 group improves the model fit slightly


#Test fit

overall.dharma <- simulateResiduals(fittedModel = overall.glmm.MRmeasured, plot = F)

plot(overall.dharma)

testDispersion(overall.glmm.MRmeasured, alternative = c("two.sided", "greater", "less"),
               plot = T, type = "DHARMa")

tab_model(overall.glmm.MRmeasured,
          show.est = TRUE) #visualizing model outputs


#### Control vs metabolic animals ----

#Create new df combining the collection data with the feeding data
controlvsmetabolic.modelling <- merge(feeding.master.May2024, PO.collection.May2024, by = "star.ID")

controlvsmetabolic.modelling <- controlvsmetabolic.modelling %>% 
  drop_na(total.mussels.consumed) %>% 
  group_by(star.ID, disc.diameter_z, trial.day, trial.block, group,
           Treatment, MR.measured.overall, lactic.acid.sampled.overall) %>% 
  summarise(total.number.mussels = sum(total.mussels.consumed))

##### MR measurements ----

#Make the glmm
controlvsmetabolic.glmm <- glmmTMB(total.number.mussels ~ MR.measured.overall +
                                     Treatment*trial.block
                                   + disc.diameter_z 
                                   + lactic.acid.sampled.overall
                                   + (1|trial.day)
                                   + (1|group),
                                   data = controlvsmetabolic.modelling, 
                                   family = poisson)

summary(controlvsmetabolic.glmm) 
#--Seems like MR measurements had a slight negative effect on mussel consumption

controlvsmetabolic.dharma <- simulateResiduals(fittedModel = controlvsmetabolic.glmm, plot = F)

plot(controlvsmetabolic.dharma) #No problems detected

testDispersion(controlvsmetabolic.glmm, alternative = c("two.sided", "greater", "less"),
               plot = T, type = "DHARMa")



#### 10-day Blocks ----

#########################i
# DAY 1 to 10
###### Adjustment ----
#########################i

adjustment.day1to10.modelling <- feeding.master.May2024 %>% 
  filter(trial.block == "Adjustment") %>% 
  drop_na(total.mussels.consumed) %>% 
  group_by(star.ID, disc.diameter_z, Treatment, water.treatment, air.treatment, group, MR.measured) %>% 
  summarise(total.number.mussels = sum(total.mussels.consumed))

#Making the GLMM
adjustment.day1to10.glmm <- glmmTMB(total.number.mussels ~ water.treatment*air.treatment +
                                      + disc.diameter_z
                                    + (1|group)
                                    + (1|MR.measured),
                                    data = adjustment.day1to10.modelling, 
                                    family = poisson)

summary(adjustment.day1to10.glmm) 
#-----Body size does not have a significant effect on mussel consumption for the adjustment period
#-----***30C air sign. diff. from intercept (higher feeding than other treatments)

adjustment.day1to10.dharma <- simulateResiduals(fittedModel = adjustment.day1to10.glmm, plot = F)

plot(adjustment.day1to10.dharma)

testDispersion(adjustment.day1to10.dharma, alternative = c("two.sided", "greater", "less"),
               plot = T, type = "DHARMa")

tab_model(adjustment.day1to10.glmm,
          show.est = TRUE)

#Create an object with the estimates and CI 
adjustment.day1to10.estimates <- plot_model(adjustment.day1to10.glmm, type = "int")


#####################################i
# Day 11 - 20
######## First half of heat stress period -----
#####################################i

heatstress.day11to20.modelling <- feeding.master.May2024 %>% 
  filter(trial.block == "First half heat stress") %>% 
  drop_na(total.mussels.consumed) %>% 
  group_by(star.ID, disc.diameter_z, water.treatment, air.treatment, group, MR.measured) %>% 
  summarise(total.number.mussels = sum(total.mussels.consumed))

#Making the GLMM
heatstress.day11to20.glmm <- glmmTMB(total.number.mussels ~ water.treatment*air.treatment +
                                       + disc.diameter_z
                                     + (1|group)
                                     + (1|MR.measured),
                                     data = heatstress.day11to20.modelling, 
                                     family = poisson)

summary(heatstress.day11to20.glmm) 
#-----No effect of body size on feeding
#-----***sign. diff. in 30C air treatment feeding

heatstress.day11to20.dharma <- simulateResiduals(fittedModel = heatstress.day11to20.glmm, plot = F)

plot(heatstress.day11to20.dharma)

testDispersion(heatstress.day11to20.dharma, alternative = c("two.sided", "greater", "less"),
               plot = T, type = "DHARMa")

tab_model(heatstress.day11to20.glmm,
          show.est = TRUE)

#Create an object with the estimates and CI 
heatstress.day11to20.estimates <- plot_model(heatstress.day11to20.glmm, type = "int")


#####################################i
# Day 21 - 30
###### Second half heat stress ----
#####################################i

heatstress.day21to30.modelling <- feeding.master.May2024 %>% 
  filter(trial.block == "Second half heat stress") %>% 
  drop_na(total.mussels.consumed) %>% 
  group_by(star.ID, disc.diameter_z, water.treatment, air.treatment, group, MR.measured) %>% 
  summarise(total.number.mussels = sum(total.mussels.consumed))

#Making the GLMM
heatstress.day21to30.glmm <- glmmTMB(total.number.mussels ~ water.treatment*air.treatment +
                                       + disc.diameter_z
                                     + (1|group)
                                     + (1|MR.measured),
                                     data = heatstress.day21to30.modelling, 
                                     family = poisson)

summary(heatstress.day21to30.glmm) 
#----We now have an effect of body size on feeding...figure out what to do about this
#----***Sign. diff. between 20C water and 15C water
#----***Sign. diff between 30C air and 25C air
#----***Sign. interaction term (water x air)

heatstress.day21to30.dharma <- simulateResiduals(fittedModel = heatstress.day21to30.glmm, plot = F)

plot(heatstress.day21to30.dharma)

testDispersion(heatstress.day21to30.dharma, alternative = c("two.sided", "greater", "less"),
               plot = T, type = "DHARMa")

tab_model(heatstress.day21to30.glmm,
          show.est = TRUE)

#Create an object with the estimates and CI 
heatstress.day21to30.estimates <- plot_model(heatstress.day21to30.glmm, type = "int")


#########################i
# Day 31 - 39
###### Recovery ----
#########################i

recovery.modelling <- feeding.master.May2024 %>% 
  filter(trial.block == "Recovery") %>% 
  drop_na(total.mussels.consumed) %>% 
  group_by(star.ID, disc.diameter_z, water.treatment, air.treatment, group, MR.measured) %>% 
  summarise(total.number.mussels = sum(total.mussels.consumed))

#Making the GLMM
recovery.glmm <- glmmTMB(total.number.mussels ~ water.treatment*air.treatment +
                           + disc.diameter_z
                         + (1|group)
                         + (1|MR.measured),
                         data = recovery.modelling, 
                         family = poisson)

summary(recovery.glmm) 
#----Body size still has significant effect...
#----***Sign. diff. between 30C air and 25C air
#----***Sign. interaction between water and air temperatures

recovery.dharma <- simulateResiduals(fittedModel = recovery.glmm, plot = F)

plot(recovery.dharma)

testDispersion(recovery.glmm, alternative = c("two.sided", "greater", "less"),
               plot = T, type = "DHARMa")

tab_model(recovery.glmm,
          show.est = TRUE)

#Create an object with the estimates and CI 
recovery.estimates <- plot_model(recovery.glmm, type = "int")


## PLOTTING ----

##Set up the data for plotting
feeding.master.plotting.2024 <- feeding.master.May2024 %>% 
  group_by(star.ID, disc.diameter_z, 
           trial.period, trial.block,
           water.treatment, air.treatment, Treatment, group) %>%
  drop_na(total.mussels.consumed) %>% 
  summarise(total.number.mussels = sum(total.mussels.consumed)) 


#Load the model estimate data
PO_modelestimates_byBlock.2024 <- read.csv("2024 data/PO_Feeding-rate_glmm-model-estimates_May2024.csv")

#Make treatments factors
PO_modelestimates_byBlock.2024$water.treatment <- as.factor(PO_modelestimates_byBlock.2024$water.treatment)
PO_modelestimates_byBlock.2024$air.treatment <- as.factor(PO_modelestimates_byBlock.2024$air.treatment)

#Add treatment column
PO_modelestimates_byBlock.2024 <- PO_modelestimates_byBlock.2024 %>% 
  mutate(Treatment = case_when(
    water.treatment =="15" & air.treatment=="25" ~ "15°C water/25°C air", 
    water.treatment =="15" & air.treatment=="30" ~ "15°C water/30°C air", 
    water.treatment =="20" & air.treatment=="25" ~ "20°C water/25°C air", 
    water.treatment =="20" & air.treatment=="30" ~ "20°C water/30°C air", 
    TRUE ~ "other"
  ))

#Make trial.block an ordered factor
feeding.master.plotting.2024$trial.block <- factor(feeding.master.plotting.2024$trial.block,
                                              ordered = TRUE,
                                              levels = c("Adjustment", "First half heat stress",
                                                         "Second half heat stress", "Recovery"))

PO_modelestimates_byBlock.2024$trial.block <- factor(PO_modelestimates_byBlock.2024$trial.block,
                                                ordered = TRUE,
                                                levels = c("Adjustment", "First half heat stress",
                                                           "Second half heat stress", "Recovery"))

### Panel B ----

#Make the baseplot with the points
Overall_PLOT_PanelB <- feeding.master.plotting.2024 %>% 
  ggplot() +
  geom_point(aes(x = Treatment, y = total.number.mussels, group = trial.block,
                 shape = trial.block, color = Treatment),
             position = position_jitterdodge(dodge.width = 1),
             alpha = 0.5) +
  scale_color_manual(name = "Treatment (°C)", 
                     values = c("#197BBD", "#05299E" , "#FC814A", "#8D0801")) +
  scale_shape_manual(name = "Trial block",
                     values = c(21, 17, 25, 8)) +
  LW_theme +
  theme(legend.position = "top",
        legend.box = "vertical",
        legend.spacing.y = unit(-0.1, "cm")) +
  labs(x = "Treatment", y = "Total # of mussels consumed per individual") +
  scale_y_continuous(breaks = c(0,1,2,3,4,5,6,7,8), limits = c(-0.6,8))

##ADDING THE ESTIMATES

Overall_PLOT_PanelB <- Overall_PLOT_PanelB +
  geom_point(data = PO_modelestimates_byBlock.2024,
             aes(x = Treatment, y = estimate, group = trial.block,
                 color = Treatment),
             position = position_dodge(width = 1), size = 4) +
  scale_color_manual(name = "Treatment (°C)", 
                     values = c("#197BBD", "#05299E" , "#FC814A", "#8D0801")) +
  #Add confidence interval for the model coefficient
  geom_errorbar(data = PO_modelestimates_byBlock.2024,
                aes(x = Treatment, y = estimate, group = trial.block,
                    ymin = conf.low, ymax = conf.high, color = Treatment),
                position = position_dodge(width = 1), width = 0.5, linewidth = 1) +
  #Add adjustment period shading
  annotate("rect", fill = "grey", alpha = 0.3, 
           xmin = 0.5, xmax = 0.75,
           ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "grey", alpha = 0.3, 
           xmin = 1.51, xmax = 1.75,
           ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "grey", alpha = 0.3, 
           xmin = 2.51, xmax = 2.75,
           ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "grey", alpha = 0.3, 
           xmin = 3.51, xmax = 3.75,
           ymin = -Inf, ymax = Inf) +
  ##Add recovery period shading
  annotate("rect", fill = "grey", alpha = 0.3, 
           xmin = 1.25, xmax = 1.49,
           ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "grey", alpha = 0.3, 
           xmin = 2.25, xmax = 2.49,
           ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "grey", alpha = 0.3, 
           xmin = 3.25, xmax = 3.49,
           ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "grey", alpha = 0.3, 
           xmin = 4.25, xmax = 4.49,
           ymin = -Inf, ymax = Inf) +
  #Add sample sizes
  geom_text(label = "n=18", x = 1, y = -0.7, size = 4) +
  geom_text(label = "n=18", x = 2, y = -0.7, size = 4) +
  geom_text(label = "n=18", x = 3, y = -0.7, size = 4) +
  geom_text(label = "n=18", x = 4, y = -0.7, size = 4)



Overall_PLOT_PanelB



###### 2023 vs 2024 -----
MS_2023vs2024 <- Overall_PLOT_PanelA + Overall_PLOT_PanelB +
  plot_layout(nrow = 2) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag.position = c(0,1),
        plot.tag = element_text(hjust = -3, vjust = 4.5))

MS_2023vs2024

#png("MS Figures/Fig3.TotalMusselsConsumed_2023vs2024.png", width = 10, height = 9, units = "in", res = 600)
#MS_2023vs2024
#dev.off()











