#################################################################################################i
# Title: Investigating the impact of contrasting water and air temperatures on the ability of
#        juvenile Pisaster ochraceus to remain attached to the aquaria (proxy for thermal stress)
#
# Author: Lydia Walton
# Last updated: 23-July-2024
################################################################################################i

#dev.off() #Clear open plots
rm(list = ls(all=T)) #Clear environment

#Set working directory
library(here)
#Packages for data manipulation
library(dplyr)
library(tidyr)
#Packages for plotting
library(ggplot2)
#Packages for statistics and modelling
library(glmmTMB)
library(DHARMa)
library(sjPlot)

# Setup ----

#------------------------------------------Load data files

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

## Data manipulation ----

#Subset the data by the variables we want 
PO.collection.May2024.sub <- PO.collection.May2024 %>% #Collection data
  subset(select = c(star.ID, ind.number, water.treatment, air.treatment, disc.diameter,
                    collection.weight))


PO.adjustment.May2024.sub <- PO.adjustment.May2024 %>% #Adjustment period (day 1 - 10)
  subset(select=c(star.ID, water.treatment, air.treatment, date, trial.period, trial.day,
                  ST.end.day,total.mussels.consumed, cumulative.mussels, thermal.stress,
                  MR.measured, lactic.acid.sampled))

PO.heat.stress.May2024.sub <- PO.heat.stress.May2024 %>% #Heat stress period (Day 11 - 29)
  subset(select=c(star.ID, water.treatment, air.treatment, date, trial.period, trial.day,
                  ST.end.day, total.mussels.consumed, cumulative.mussels, thermal.stress,
                  MR.measured, lactic.acid.sampled))

PO.recovery.May2024.sub <- PO.recovery.May2024 %>% #Heat stress period (Day 11 - 29)
  subset(select=c(star.ID, water.treatment, air.treatment, date, trial.period, trial.day,
                  ST.end.day,total.mussels.consumed, cumulative.mussels, thermal.stress,
                  MR.measured, lactic.acid.sampled))

#Merge dfs together
attachment.master.May2024 <-rbind(PO.adjustment.May2024.sub, PO.heat.stress.May2024.sub)

attachment.master.May2024 <-rbind(attachment.master.May2024, PO.recovery.May2024.sub)

attachment.master.May2024 <- merge(PO.collection.May2024.sub, attachment.master.May2024,
                                by = c("star.ID", "water.treatment", "air.treatment"))

#Make sure all of the variables are the right structure
attachment.master.May2024$water.treatment <- as.factor(attachment.master.May2024$water.treatment)
attachment.master.May2024$air.treatment <- as.factor(attachment.master.May2024$air.treatment)

##RESCALE body size for modelling
attachment.master.May2024$disc.diameter_z <- arm::rescale(attachment.master.May2024$disc.diameter)

attachment.master.May2024 <- attachment.master.May2024 %>% 
  mutate(Treatment = case_when(
    water.treatment =="15" & air.treatment=="25" ~ "15°C water/25°C air", 
    water.treatment =="15" & air.treatment=="30" ~ "15°C water/30°C air", 
    water.treatment =="20" & air.treatment=="25" ~ "20°C water/25°C air", 
    water.treatment =="20" & air.treatment=="30" ~ "20°C water/30°C air", 
    TRUE ~ "other"
  ))

#New column for trial "blocks"
attachment.master.May2024 <- attachment.master.May2024 %>% 
  mutate(trial.block = case_when(
    trial.day < 11 ~ "Adjustment", 
    trial.day > 10 & trial.day < 21 ~ "First half heat stress", 
    trial.day > 20 & trial.day < 31 ~ "Second half heat stress",
    trial.day > 30 ~ "Recovery",
    TRUE ~ "OTHER"
  ))

#Add group column based on odd or even ind.number
attachment.master.May2024 <- attachment.master.May2024 %>% 
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
  drop_na(thermal.stress)


## Thermal stress ----

### Statistics ----

# We want to test the difference in attachment (yes or no) across treatments

#GLMM with a binomial distribution

attachment.glmm.MRmeasured <- glmmTMB(thermal.stress ~ trial.block*Treatment +
                             disc.diameter_z + MR.measured +
                             (1|star.ID) +
                             (1|group),
                           data = attachment.master.May2024, family = binomial)

summary(attachment.glmm.MRmeasured)


#MR handling as random effect instead
attachment.glmm <- glmmTMB(thermal.stress ~ trial.block*Treatment +
                             disc.diameter_z +
                             (1|MR.measured) +
                             (1|star.ID) +
                             (1|group),
                           data = attachment.master.May2024, family = binomial)

summary(attachment.glmm)

#Nested random effects of individuals within MO2 group
attachment.glmm.nested <- glmmTMB(thermal.stress ~ trial.block*Treatment +
                             disc.diameter_z +
                             (1|MR.measured/star.ID),
                           data = attachment.master.May2024, family = binomial)

summary(attachment.glmm.nested)

#Compare model fit
AIC(attachment.glmm.MRmeasured,attachment.glmm, attachment.glmm.nested)
##Including MR.measured gives a slightly better fit - there seems to be a slight effect of handling

#-----Recovery had sign. higher mussel consumption compared to adjustment period
attachment.dharma <- simulateResiduals(fittedModel = attachment.glmm.nested, plot = F)

plot(attachment.dharma) #No problems detected

testDispersion(attachment.glmm.nested, alternative = c("two.sided", "greater", "less"),
               plot = T, type = "DHARMa")

tab_model(attachment.glmm.nested,
          show.est = TRUE) #visualizing model outputs


### Average # detached stars ----
ave.detached.stars.pertreatment <- attachment.master.May2024 %>% 
  filter(thermal.stress == "1") %>% 
  group_by(Treatment, trial.period, trial.block) %>% 
  summarise(number.detached = length(thermal.stress), percent.detached = number.detached/180*100) 

ave.attached.stars.pertreatment <- attachment.master.May2024  %>% 
  filter(thermal.stress == "0") %>% 
  group_by(Treatment, trial.period, trial.block) %>% 
  summarise(number.attached = length(thermal.stress))

### Plots ----

#Plot the data to see what it looks like

#Order the trial blocks
attachment.master.May2024$trial.block <- factor(attachment.master.May2024$trial.block,
                                                ordered = TRUE,
                                                levels = c("Adjustment", "First half heat stress",
                                                           "Second half heat stress", "Recovery"))

#Filter the data for stars that were recorded as being thermally stressed
attachment.plotting <- attachment.master.May2024 %>% 
  mutate(attached = case_when(
    thermal.stress == "0" ~ "Atttached",
    thermal.stress == "1" ~ "Not attached",
    TRUE ~ "other")) %>% 
 # filter(attached == "Not attached") %>% 
  filter(trial.day > 12)

#Make df for panel labels (A-D)
label_text <- data.frame(
  label = c("A", "B", "C", "D"),
  Treatment = c("15°C water/25°C air", 
            "15°C water/30°C air", 
            "20°C water/25°C air", 
            "20°C water/30°C air"))

#Make the plot
attachment.barplot <- attachment.plotting %>% 
  ggplot() +
  geom_bar(aes(x = trial.day, fill = attached)) +
  scale_fill_manual(name = "Attachment to aquaria",
                    values = c("#9AD5CA" ,"#90323D")) +
  LW_theme +
  facet_wrap(~Treatment) +
  scale_x_continuous(breaks = seq(13, 38, by = 2), expand = c(0.01,0.01)) +
  scale_y_continuous(breaks = seq(1,18, by = 2), expand = c(0.01, 0.5)) +
  labs(x = "Day of experiment", y = "Number of animals") +
  theme(legend.position = "bottom") +
  #Add shading for recovery period
  annotate("rect", fill = "grey", alpha = 0.65, 
           xmin = 30.5, xmax = 38.5,
           ymin = -Inf, ymax = Inf) +
  #Add panel labels
  geom_text(data = label_text,
            mapping = aes(x = -Inf, y = Inf, label = label, size = 4),
            hjust   = -1,
            vjust   = 2,
            show.legend = FALSE)

attachment.barplot

#png("MS Figures/Fig4.Attachment-to-aquaria_May2024.png", width = 10, height = 9, units = "in", res = 600)
#attachment.barplot
#dev.off()




