#################################################################################################i
# Title: Investigating the impact of contrasting water and air temperatures on the metabolic rate
#        (estimated through aquatic oxygen consumption) of juvenile Pisaster ochraceus
#
# Author: Lydia Walton
# Last updated: 13-Dec-2024
################################################################################################i

#dev.off() #Clear open plots
rm(list = ls(all=T)) #Clear environment

#Set working directory
library(here)
#Packages for data manipulation
library(dplyr)
library(tidyr)
#Data visualization
library(ggplot2)
library(ggpmisc)
library(tibble)
library(kableExtra)
library(sjPlot)
library(RColorBrewer)
library(extrafont)
library(ggpubr)
library(patchwork)
#Modeling
library(quantreg)
library(AICcmodavg)
library(glmmTMB)
library(EnvStats)
library(broom)
library(DHARMa)
library(emmeans)
library(rstatix)
library(nlme)


# Setup ----

#------------------------------------------Load data files

PO.collection.May2024 <- read.csv("2024 data/PO_Collection_May-2024.csv")
PO.adjustment.May2024 <- read.csv("2024 data/PO_Adjustment_May-2024.csv")
PO.heat.stress.May2024 <- read.csv("2024 data/PO_Heatstress_May-2024.csv")
PO.recovery.May2024 <- read.csv("2024 data/PO_Recovery_May-2024.csv")

#Metabolic data
M1.metadata <- read.csv("2024 data/PO_M1-metadata_21-May-24.csv")

M1.O2.values <- read.csv("2024 data/PO_M1-O2-values_21-May-24.csv")

M2.metadata <- read.csv("2024 data/PO_M2-metadata_31-May-24.csv")

M2.O2.values <- read.csv("2024 data/PO_M2-O2-values_31-May-24.csv")

M3.metadata <- read.csv("2024 data/PO_M3-metadata_09-June-24.csv")

M3.O2.values <- read.csv("2024 data/PO_M3-O2-values_09-June-24.csv")

M4.metadata <- read.csv("2024 data/PO_M4-metadata_19-June-24.csv")

M4.O2.values <- read.csv("2024 data/PO_M4-O2-values_19-June-24.csv")


###Create theme for plotting

#Lydia's theme
LW_theme <- theme_classic() +
  theme(axis.title.x = element_text(vjust = -1, size = 12),
        axis.title.y = element_text(vjust = 2, size = 11),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 10),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 10),
        panel.background = element_rect(fill = NA, color = "black"))


## Data manipulation ----

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

#Add column for Treatment
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

#-------------------------Subset metabolic spreadsheets for merging
#M1 metadata
M1.metadata.sub <- M1.metadata %>% 
  subset(select = c(water.treatment, air.treatment, star.ID, ind.number, measurement.number, 
                    vol.water.L, wet.weight.g)) %>% 
  drop_na(ind.number)

#M2 metadata
M2.metadata.sub <- M2.metadata %>% 
  subset(select = c(water.treatment, air.treatment, star.ID, ind.number, measurement.number, 
                    vol.water.L, wet.weight.g)) %>% 
  drop_na(ind.number)

#M3 metadata
M3.metadata.sub <- M3.metadata %>% 
  subset(select = c(water.treatment, air.treatment, star.ID, ind.number, measurement.number, 
                    vol.water.L, wet.weight.g)) %>% 
  drop_na(ind.number)

#M4 metadata
M4.metadata.sub <- M4.metadata %>% 
  subset(select = c(water.treatment, air.treatment, star.ID, ind.number, measurement.number, 
                    vol.water.L, wet.weight.g)) %>% 
  drop_na(ind.number)

#----------------O2 value datasheets
#M1 O2 values
M1.O2.values.sub <- M1.O2.values %>% 
  subset(select = c(water.treatment, air.treatment, star.ID, ind.number, MO2.type, round, channel, 
                    salinity.ppm, unconverted.rate, AFDM.MR.mLO2.hr.g, 
                    rsqrt.MR.weight, percent.drop.O2)) %>% 
  drop_na(ind.number)

#M2 O2 values
M2.O2.values.sub <- M2.O2.values %>% 
  subset(select = c(water.treatment, air.treatment, star.ID, ind.number, MO2.type, round, channel, 
                    salinity.ppm, unconverted.rate, AFDM.MR.mLO2.hr.g, 
                    rsqrt.MR.weight, percent.drop.O2)) %>% 
  drop_na(ind.number)

#M3 O2 values
M3.O2.values.sub <- M3.O2.values %>% 
  subset(select = c(water.treatment, air.treatment, star.ID, ind.number, MO2.type, round, channel, 
                    salinity.ppm, unconverted.rate, AFDM.MR.mLO2.hr.g, 
                    rsqrt.MR.weight, percent.drop.O2)) %>% 
  drop_na(ind.number)

#M4 O2 values
M4.O2.values.sub <- M4.O2.values %>% 
  subset(select = c(water.treatment, air.treatment, star.ID, ind.number, MO2.type, round, channel, 
                    salinity.ppm, unconverted.rate, AFDM.MR.mLO2.hr.g, 
                    rsqrt.MR.weight, percent.drop.O2)) %>% 
  drop_na(ind.number)


#---------------------------Merge metabolic dataframes into one master
#M1
M1.data <- merge(M1.metadata.sub, M1.O2.values.sub, by = c("water.treatment", "air.treatment",
                                                           "star.ID",
                                                           "ind.number"))
#M2
M2.data <- merge(M2.metadata.sub, M2.O2.values.sub, by = c("water.treatment", "air.treatment",
                                                           "star.ID",
                                                           "ind.number"))
#M3
M3.data <- merge(M3.metadata.sub, M3.O2.values.sub, by = c("water.treatment", "air.treatment",
                                                           "star.ID",
                                                           "ind.number"))

#M4
M4.data <- merge(M4.metadata.sub, M4.O2.values.sub, by = c("water.treatment", "air.treatment",
                                                           "star.ID",
                                                           "ind.number"))

#Bind the metabolic measurements together

Metabolic.master <- rbind(M1.data, M2.data)

Metabolic.master <-rbind(Metabolic.master, M3.data)

Metabolic.master.aq.aer <-rbind(Metabolic.master, M4.data)


#---------------------------Cleaning columns
#As factors
Metabolic.master.aq.aer$air.treatment <- as.factor(Metabolic.master.aq.aer$air.treatment)
Metabolic.master.aq.aer$water.treatment <- as.factor(Metabolic.master.aq.aer$water.treatment)

#Change the MO2 data to absolute values (no negative signs)
Metabolic.master.aq.aer$AFDM.MR.mLO2.hr.g <- abs(Metabolic.master.aq.aer$AFDM.MR.mLO2.hr.g)

#Filter out everything except First.Aquatic measurements
Metabolic.master <- Metabolic.master.aq.aer %>% 
  filter(MO2.type == "First Aquatic")

#Make a new column "Treatment" for plotting
Metabolic.master <- Metabolic.master %>% 
  mutate(Treatment = case_when(
    water.treatment =="15" & air.treatment=="25" ~ "15°C water/25°C air", 
    water.treatment =="15" & air.treatment=="30" ~ "15°C water/30°C air",  
    water.treatment =="20" & air.treatment=="25" ~ "20°C water/25°C air", 
    water.treatment =="20" & air.treatment=="30" ~ "20°C water/30°C air",  
    TRUE ~ "other"
  ))

Metabolic.master <- Metabolic.master %>% 
  mutate(trial.block = case_when(
    measurement.number == "M1" ~ "Adjustment", 
    measurement.number == "M2" ~ "First half heat stress", 
    measurement.number == "M3" ~ "Second half heat stress",
    measurement.number == "M4" ~ "Recovery",
    TRUE ~ "OTHER"
  ))

#Add new labels for the facet wrap Water temp and Air temp
water.labels <- c("15°C water treatment", "20°C water treatment")
names(water.labels) <- c("15", "20")

air.labels <- c("25°C air treatment", "30°C air treatment")
names(air.labels) <- c("25", "30")

measurement.labels <- c("Adjustment", "First half heat stress", 
                        "Second half heat stress", "Recovery")
names(measurement.labels) <- c("M1", "M2", "M3", "M4")


# Aquatic MO2 ----

## Statistics ----

#New df for modelling (add feeding and thermal stress data to metabolic data)

##Summarise PO.master
metabolic.modelling <- feeding.master.May2024 %>%
  drop_na(total.mussels.consumed) %>% 
  group_by(water.treatment, air.treatment, star.ID, ind.number, disc.diameter_z,
           trial.period, trial.block, MR.measured, lactic.acid.sampled) %>% 
  summarise(total.mussels = sum(total.mussels.consumed))

#merge with metabolic.master
metabolic.modelling <- merge(metabolic.modelling, Metabolic.master, 
                             by = c("water.treatment", "air.treatment", "star.ID", "ind.number",
                                    "trial.block"))

#Create another df for plotting (will use later)
metabolic.plotting <- metabolic.modelling


### Linear mixed model ----
#Get rid of MR with an R-squared less than 0.98
metabolic.modelling.rsqrt98 <- metabolic.modelling %>% 
  filter(rsqrt.MR.weight > 0.98) %>% 
  drop_na(AFDM.MR.mLO2.hr.g)

### modelling changes in metabolic rate
#channel isn't needed, if you look in the summary table including channel explains zero variance in the data, which means it is equivalent to random
#weights models the variance structure of each treatment, and may improve the model fit
metabolic.modelling.rsqrt98$Treatment<-as.factor(metabolic.modelling.rsqrt98$Treatment)
metabolic.modelling.rsqrt98$Treatment<-relevel(metabolic.modelling.rsqrt98$Treatment,
                                               ref="15°C water/30°C air")

metabolic.modelling.rsqrt98$MR.measured <-as.factor(metabolic.modelling.rsqrt98$MR.measured)


MR.lme <- lme(AFDM.MR.mLO2.hr.g ~ Treatment * measurement.number + disc.diameter_z
              + total.mussels,
              random = (~1|star.ID), weights=varIdent(form = ~1|Treatment),
              data = metabolic.modelling.rsqrt98)


MR.lme.no.weights <- lme(AFDM.MR.mLO2.hr.g ~ Treatment * measurement.number +
                           disc.diameter_z + total.mussels, 
                         random = (~1|star.ID),
                         data = metabolic.modelling.rsqrt98)

#weights does improve the model fit
AIC(MR.lme,MR.lme.no.weights)

summary(MR.lme)
qqnorm(MR.lme)
plot(MR.lme)

tab_model(MR.lme,
          show.est = TRUE)

library(lsmeans)

#Look at differences within each measurement
variables<-lsmeans(MR.lme, ~Treatment|measurement.number)
pairs(variables,adjust = "tukey")

#Look at differences within each treatment
variables.2<-lsmeans(MR.lme, ~measurement.number|Treatment)
pairs(variables.2,adjust = "tukey")


## Plotting ----

#Get rid of MR with an R-squared less than 0.98
metabolic.plotting.rsqrt98 <- metabolic.plotting %>% 
  filter(rsqrt.MR.weight > 0.98) %>% 
  drop_na(AFDM.MR.mLO2.hr.g)


#Sample sizes for each measurement
MR.sample.size.rsqrt98 <- metabolic.modelling.rsqrt98 %>% 
  group_by(Treatment, trial.block, measurement.number) %>% 
  summarise(sample.size = length(AFDM.MR.mLO2.hr.g))

MR.sample.size.rsqrt98

#Make df for panel labels (A-D)
label_text <- data.frame(
  label = c("A", "B", "C", "D"),
  measurement.number = c("M1", "M2", "M3", "M4"))

####PANEL A - Adjustment ----
metabolic.plotting.rsqrt98.adjustment <- metabolic.plotting.rsqrt98 %>% 
  filter(measurement.number == "M1")


#----------Boxplot
MR.boxplot.rsqrt98.adj <- metabolic.plotting.rsqrt98.adjustment %>% 
  ggplot() +
  geom_boxplot(aes(x = Treatment, y = AFDM.MR.mLO2.hr.g, fill = Treatment), outlier.shape = NA) +
  scale_fill_manual(values = c("#197BBD", "#05299E" , "#FC814A", "#8D0801")) +
  facet_wrap(~measurement.number,
             labeller = labeller(measurement.number = measurement.labels)) +
  LW_theme +
  theme(legend.position = "top",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(x = "Treatment", y = expression(dot(M)*O[2] ~ "per g of AFDM")) +
  ylim(0.1, 0.45) +
  stat_n_text(mapping = (aes(x = Treatment, y = AFDM.MR.mLO2.hr.g)),
    y.pos = 0.1, color = "black") +
  #Add panel label
  geom_text(label= "A", x = 0.6, y = 0.45, 
            size = 5, color = "black")+
  #Add panel shading
  annotate("rect", fill = "grey", alpha = 0.3, 
         xmin = -Inf, xmax = Inf,
         ymin = -Inf, ymax = Inf) +
  #Add Tukey line
  geom_segment(aes(x = 1, y = 0.41, xend = 2, yend = 0.41)) +
  #Add sig.level
  geom_text(label = "*", x = 1.5, y = 0.415, size = 5)

MR.boxplot.rsqrt98.adj

####PANEL B - Heat stress ----
metabolic.plotting.rsqrt98.firstheat <- metabolic.plotting.rsqrt98 %>% 
  filter(measurement.number == "M2")


#----------Boxplot
MR.boxplot.rsqrt98.heat1 <- metabolic.plotting.rsqrt98.firstheat %>% 
  ggplot() +
  geom_boxplot(aes(x = Treatment, y = AFDM.MR.mLO2.hr.g, fill = Treatment), outlier.shape = NA) +
  scale_fill_manual(values = c("#197BBD", "#05299E" , "#FC814A", "#8D0801")) +
  facet_wrap(~measurement.number,
             labeller = labeller(measurement.number = measurement.labels)) +
  LW_theme +
  theme(legend.position = "top",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  labs(x = "Treatment", y = expression(dot(M)*O[2] ~ "per g of AFDM")) +
  ylim(0.1, 0.45) +
  stat_n_text(mapping = (aes(x = Treatment, y = AFDM.MR.mLO2.hr.g)),
              y.pos = 0.1, color = "black") +
  #Add panel label
  geom_text(label= "B", x = 0.6, y = 0.45, 
            size = 5, color = "black") +
  #Add Tukey line
  geom_segment(aes(x = 2, y = 0.43, xend = 4, yend = 0.43)) +
  geom_segment(aes(x = 1, y = 0.45, xend = 4, yend = 0.45)) +
  geom_segment(aes(x = 1, y = 0.35, xend = 3, yend = 0.35)) +
  #Add sig.level
  geom_text(label = "*", x = 3, y = 0.435, size = 5) +
  geom_text(label = "***", x = 2.5, y = 0.455, size = 5) +
  geom_text(label = "***", x = 2, y = 0.355, size = 5)

MR.boxplot.rsqrt98.heat1


####PANEL C - Heat stress ----
metabolic.plotting.rsqrt98.secondheat <- metabolic.plotting.rsqrt98 %>% 
  filter(measurement.number == "M3")


#----------Boxplot
MR.boxplot.rsqrt98.heat2 <- metabolic.plotting.rsqrt98.secondheat %>% 
  ggplot() +
  geom_boxplot(aes(x = Treatment, y = AFDM.MR.mLO2.hr.g, fill = Treatment), outlier.shape = NA) +
  scale_fill_manual(values = c("#197BBD", "#05299E" , "#FC814A", "#8D0801")) +
  facet_wrap(~measurement.number,
             labeller = labeller(measurement.number = measurement.labels)) +
  LW_theme +
  theme(legend.position = "top") +
  labs(x = "Treatment", y = expression(dot(M)*O[2] ~ "per g of AFDM")) +
  ylim(0.1, 0.45) +
  stat_n_text(mapping = (aes(x = Treatment, y = AFDM.MR.mLO2.hr.g)),
              y.pos = 0.1, color = "black", family = "A") +
  #Add panel label
  geom_text(label= "C", x = 0.6, y = 0.45, 
            size = 5, family = "A", color = "black") +
  #Add Tukey line
  geom_segment(aes(x = 2, y = 0.37, xend = 3, yend = 0.37)) +
  geom_segment(aes(x = 1, y = 0.4, xend = 3, yend = 0.4)) +
  #Add sig.level
  geom_text(label = "*", x = 2.5, y = 0.375, size = 5) +
  geom_text(label = "***", x = 2, y = 0.405, size = 5)

MR.boxplot.rsqrt98.heat2

####PANEL D - Recovery ----
metabolic.plotting.rsqrt98.rec <- metabolic.plotting.rsqrt98 %>% 
  filter(measurement.number == "M4")


#----------Boxplot
MR.boxplot.rsqrt98.rec <- metabolic.plotting.rsqrt98.rec %>% 
  ggplot() +
  geom_boxplot(aes(x = Treatment, y = AFDM.MR.mLO2.hr.g, fill = Treatment), outlier.shape = NA) +
  scale_fill_manual(values = c("#197BBD", "#05299E" , "#FC814A", "#8D0801")) +
  facet_wrap(~measurement.number,
             labeller = labeller(measurement.number = measurement.labels)) +
  LW_theme +
  theme(legend.position = "top",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  labs(x = "Treatment", y = expression(dot(M)*O[2] ~ "per g of AFDM")) +
  ylim(0.1, 0.45) +
  stat_n_text(mapping = (aes(x = Treatment, y = AFDM.MR.mLO2.hr.g)),
              y.pos = 0.1, color = "black", family = "A") +
  #Add panel label
  geom_text(label= "D", x = 0.6, y = 0.45, 
            size = 5, family = "A", color = "black") +
  #Add panel shading
  annotate("rect", fill = "grey", alpha = 0.3, 
           xmin = -Inf, xmax = Inf,
           ymin = -Inf, ymax = Inf)  

MR.boxplot.rsqrt98.rec

#Plot for MS ----
MS.metabolicrates.plot <- MR.boxplot.rsqrt98.adj + MR.boxplot.rsqrt98.heat1 + MR.boxplot.rsqrt98.heat2 + MR.boxplot.rsqrt98.rec +
  plot_layout(nrow = 2, ncol = 2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top") 

MS.metabolicrates.plot

#png("MS Figures/Fig5.Aquatic-MO2-Boxplot.png", width = 10, height = 9, units = "in", res = 600)
#MS.metabolicrates.plot
#dev.off()


# Aerial MO2 ----

metabolic.master.aerial <- Metabolic.master.aq.aer %>% 
  filter(MO2.type == "Aerial") %>% 
  filter(measurement.number != "M3") %>% 
  mutate(Treatment = case_when(
    water.treatment =="15" & air.treatment=="25" ~ "15°C water/25°C air", 
    water.treatment =="15" & air.treatment=="30" ~ "15°C water/30°C air",  
    water.treatment =="20" & air.treatment=="25" ~ "20°C water/25°C air", 
    water.treatment =="20" & air.treatment=="30" ~ "20°C water/30°C air",  
    TRUE ~ "other"
  )) %>% 
  mutate(trial.block = case_when(
    measurement.number == "M1" ~ "Adjustment", 
    measurement.number == "M2" ~ "First half heat stress", 
    measurement.number == "M3" ~ "Second half heat stress",
    measurement.number == "M4" ~ "Recovery",
    TRUE ~ "OTHER"
  ))

## Summary table ----
#Make table with the min, max and mean aerial O2 consumption
Aerial.MR.summarytable <- metabolic.master.aerial %>% 
  group_by(Treatment, trial.block) %>% 
  summarise(Mean = mean(percent.drop.O2),
            SD = sd(percent.drop.O2),
            Min = min(percent.drop.O2),
            Max = max(percent.drop.O2))

Aerial.MR.summarytable

## Statistics ----

#Compare % drop in O2 across treatments 

## ANOVA

aerial.aov <- aov(percent.drop.O2 ~ Treatment*measurement.number,
    data = metabolic.master.aerial)

shapiro.test(aerial.aov$residuals) #almost Normal
hist(aerial.aov$residuals)
leveneTest(percent.drop.O2 ~ Treatment*measurement.number,
           data = metabolic.master.aerial)

summary(aerial.aov)

TukeyHSD(aerial.aov)

#non-parametric version
kruskal.test(percent.drop.O2 ~ Treatment,
             data = metabolic.master.aerial)

kruskal.test(percent.drop.O2 ~ measurement.number,
             data = metabolic.master.aerial) #significant

## Plotting ----

#Make df for panel labels (A-D)
label_text <- data.frame(
  label = c("A", "B"),
  measurement.number = c("M1", "M2"))

#make variables factors for plotting
metabolic.aerial.plotting <- metabolic.master.aerial

metabolic.aerial.plotting$Treatment <- as.factor(metabolic.aerial.plotting$Treatment)

#----------Boxplot
MR.boxplot.aerial <- metabolic.master.aerial %>% 
  ggplot() +
  geom_boxplot(aes(x = Treatment, y = percent.drop.O2, fill = Treatment)) +
  scale_fill_manual(values = c("#197BBD", "#05299E" , "#FC814A", "#8D0801")) +
  facet_wrap(~measurement.number,
             labeller = labeller(measurement.number = measurement.labels)) +
  LW_theme +
  theme(legend.position = "top") +
  labs(x = "Treatment", y = "% drop in oxygen in air") +
  stat_n_text(mapping = (aes(x = Treatment, y = percent.drop.O2)),
              y.pos = 0.1, color = "black", family = "A") +
  #Add panel labels
  geom_text(data = label_text,
            mapping = aes(x = -Inf, y = Inf, label = label, family = "A", size = 4),
            hjust   = -1,
            vjust   = 2,
            show.legend = FALSE)

MR.boxplot.aerial

#png("SuppMat Figures/FigS2.Aerial-percentO2.png", width = 10, height = 9, units = "in", res = 600)
#MR.boxplot.aerial
#dev.off()









