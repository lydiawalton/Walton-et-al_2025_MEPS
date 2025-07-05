##############################################################################################
# Title: Investigating the impact of contrasting water and air temperatures on the feeding
#        activity of juvenile Pisaster ochraceus (Part 2)
#
# Author: Lydia Walton
# Last updated: 19-Mar-2024
############################################################################################

dev.off() #Clear open plots
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
library(changepoint)
library(car)
library(rstatix)
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


##Load data files
seastar.collection <- read.csv("2023 Data/Pisaster_Collections_June2023.csv")
seastar.adj.heat <- read.csv("2023 Data/Pisaster_HeatstressPeriod_June2023.csv")
seastar.recovery <- read.csv("2023 Data/Pisaster_RecoveryPeriod_June2023.csv")


###Create theme for plotting
#Font change to Times New Roman as "A"
#windowsFonts(A = windowsFont("Times New Roman"))

#Lydia's theme
LW_theme <- theme_classic() +
  theme(#text = element_text(family = "A"),
        axis.title.x = element_text(vjust = -1, size = 12),
        axis.title.y = element_text(vjust = 2, size = 11),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        panel.background = element_rect(fill = NA, color = "black"))


##----------------Data cleaning and manipulation
#Subset the data by the variables we want 
subset.collection <- seastar.collection %>% #Collection data
  subset(select = c(star.ID, arm.length, disc.diameter, collection.weight))


subset.adj.heat <- seastar.adj.heat %>% #Adj day (day9) + heat stress
  subset(select=c(star.ID, water.treatment, air.treatment, trial.period, heat.trial.day, total.mussels.consumed, cum.mussels, trial.type, group))


subset.recovery <- seastar.recovery %>% #Recovery
  subset(select=c(star.ID, water.treatment, air.treatment, trial.period, rec.trial.day, rec.mussels.consumed, cum.mussels, trial.type, group))

#Change the column names to match
colnames(subset.adj.heat) <- c("star.ID", "water.treatment", "air.treatment", "trial.period", "day.of.trial", "mussels.consumed","cumulative.mussels", "trial.type", "group")

colnames(subset.recovery) <- c("star.ID", "water.treatment", "air.treatment", "trial.period", "day.of.trial", "mussels.consumed","cumulative.mussels", "trial.type", "group")

#Merge dfs together
master.number <- rbind(subset.adj.heat, subset.recovery)

#Add the collection data
master.number <- merge(master.number, subset.collection, by = "star.ID")

#Make sure all of the variables are the right structure
master.number$water.treatment <- as.factor(master.number$water.treatment)
master.number$air.treatment <- as.factor(master.number$air.treatment)

#Remove #67 (4 armed star) and the sea stars that died later
master.number <- master.number %>% 
  filter(star.ID != "SB_2023_ind67", star.ID != "SB_2023_ind9", star.ID != "SB_2023_ind11", star.ID != "SB_2023_ind18", star.ID != "SB_2023_ind19", star.ID != "SB_2023_ind21", star.ID != "SB_2023_ind23", star.ID != "SB_2023_ind53", star.ID != "SB_2023_ind58")

#Remove NAs
master.number <- na.omit(master.number)

##RESCALE body size for modelling
master.number$disc.diameter_z <- arm::rescale(master.number$disc.diameter)


##---------------------------------------------------------
# Looking at the cumulativeulative consumption of mussels by PO

#Add new labels for the facet wrap Water temp and Air temp
water.labels <- c("15°C water treatment", "20°C water treatment")
names(water.labels) <- c("15", "20")

air.labels <- c("20°C air treatment", "25°C air treatment", "30°C air treatment")
names(air.labels) <- c("20", "25", "30")


###############################
#
# CHANGE POINT ANALYSIS
#
###############################

#Make a subset of the data that only includes the variables you need (Day of trial and mussels consumed)
total.cumulative.consumption <- master.number %>% 
  filter(trial.type == "Heat stress" | trial.type == "Recovery") %>% 
  group_by(water.treatment, air.treatment, trial.type, day.of.trial) %>% 
  summarise(total.cumulative.number.mussels = sum(cumulative.mussels),
            total.number.mussels = sum(mussels.consumed))

changepoint.master <- total.cumulative.consumption %>% 
  subset(select = c(air.treatment, water.treatment, day.of.trial, total.cumulative.number.mussels, total.number.mussels))

changepoint.master$day.of.trial <- as.numeric(changepoint.master$day.of.trial)
changepoint.master$total.cumulative.number.mussels <- as.numeric(changepoint.master$total.cumulative.number.mussels)
changepoint.master$total.cumulative.number.mussels <- as.numeric(changepoint.master$total.number.mussels)

#Subset the data by treatment group 

#-------------30 air and 15 water---------------------
cp.30air15water.subset <- changepoint.master %>% 
  filter(air.treatment == "30") %>% 
  filter(water.treatment == "15")

obj2.30A15W <- cpt.mean(cp.30air15water.subset$total.number.mussels, method = "AMOC")

plot(obj2.30A15W)

print(obj2.30A15W)

cp.30air15water.subset[17,] ##Change point detected on day 18 of the trial

#-----------30 air and 20 water------------------------
cp.30air20water.subset <- changepoint.master %>% 
  filter(air.treatment == "30") %>% 
  filter(water.treatment == "20")

obj2.30A20W <- cpt.mean(cp.30air20water.subset$total.number.mussels, method = "AMOC")

plot(obj2.30A20W)

print(obj2.30A20W)

cp.30air20water.subset[12,] ##Change point detected on day 13

#-----------25 air and 15 water------------------------
cp.25air15water.subset <- changepoint.master %>% 
  filter(air.treatment == "25") %>% 
  filter(water.treatment == "15")

obj2.25A15W <- cpt.mean(cp.25air15water.subset$total.number.mussels, method = "AMOC")

plot(obj2.25A15W)

print(obj2.25A15W)

cp.25air15water.subset[14,] ##Change point detected on day 15

#--------------25 air and 20 water-------------------------
cp.25air20water.subset <- changepoint.master %>% 
  filter(air.treatment == "25") %>% 
  filter(water.treatment == "20")

obj2.25A20W <- cpt.mean(cp.25air20water.subset$total.number.mussels, method = "AMOC")

plot(obj2.25A20W)

print(obj2.25A20W)

cp.25air20water.subset[11,] ##Change point detected on day 12

#------------20 air and 15 water--------------------------
cp.20air15water.subset <- changepoint.master %>% 
  filter(air.treatment == "20") %>% 
  filter(water.treatment == "15")

obj2.20A15W <- cpt.mean(cp.20air15water.subset$total.number.mussels, method = "AMOC")

plot(obj2.20A15W)

print(obj2.20A15W)

cp.20air15water.subset[12,] ##Change point detected on day 13

#------------20 air and 20 water-------------------------

cp.20air20water.subset <- changepoint.master %>% 
  filter(air.treatment == "20") %>% 
  filter(water.treatment == "20")

obj2.20A20W <- cpt.mean(cp.20air20water.subset$total.number.mussels, method = "AMOC")

plot(obj2.20A20W)

print(obj2.20A20W)

cp.20air20water.subset[9,] ##Change point detected on day 10

###THE SAME CHANGEPOINT DAYS WERE DETECTED USING BOTH cumulativeULATIVE AND TOTAL MUSSEL CONSUMPTION

###################################################
## PLOTTING
###################################################

#Add a new column with the combined treatment types
total.cumulative.consumption.PLOTTING <- total.cumulative.consumption %>%
  mutate(Treatment = case_when(
    water.treatment =="15" & air.treatment=="20" ~ "15°C water/20°C air",
    water.treatment =="15" & air.treatment=="25" ~ "15°C water/25°C air", 
    water.treatment =="15" & air.treatment=="30" ~ "15°C water/30°C air", 
    water.treatment =="20" & air.treatment=="20" ~ "20°C water/20°C air",
    water.treatment =="20" & air.treatment=="25" ~ "20°C water/25°C air", 
    water.treatment =="20" & air.treatment=="30" ~ "20°C water/30°C air", 
    TRUE ~ "other"
  ))

## Add new column delineating where the changepoints are for each treatment
total.cumulative.consumption.PLOTTING <- total.cumulative.consumption.PLOTTING %>%
  mutate(changepoint = case_when(
  Treatment == "15°C water/20°C air" &  day.of.trial < 13 ~ "Before",
  Treatment == "15°C water/20°C air" &  day.of.trial > 12 ~ "After",
  Treatment == "15°C water/25°C air" &  day.of.trial < 15 ~ "Before",
  Treatment == "15°C water/25°C air" &  day.of.trial > 14 ~ "After",
  Treatment == "15°C water/30°C air" &  day.of.trial < 18 ~ "Before",
  Treatment == "15°C water/30°C air" &  day.of.trial > 17 ~ "After",
  Treatment == "20°C water/30°C air" &  day.of.trial < 13 ~ "Before",
  Treatment == "20°C water/30°C air" &  day.of.trial > 12 ~ "After",
  Treatment == "20°C water/25°C air" &  day.of.trial < 12 ~ "Before",
  Treatment == "20°C water/25°C air" &  day.of.trial > 11 ~ "After",
  Treatment == "20°C water/20°C air" &  day.of.trial < 10 ~ "Before",
  Treatment == "20°C water/20°C air" &  day.of.trial > 9 ~ "After",
  TRUE ~ "other"
)) 

#Change changepoint column to factor
total.cumulative.consumption.PLOTTING$changepoint <- factor(total.cumulative.consumption.PLOTTING$changepoint,
                                            ordered = TRUE,
                                            levels = c("Before", "After"))

##Add new rows to remove the gap in the data when plotting
##15 water and 20 air
total.cumulative.consumption.PLOTTING[nrow(total.cumulative.consumption.PLOTTING) + 1,] <- 
  list("15", "20", "Heat stress", 13, 76, 3, "15°C water/20°C air", "Before")
##15 water and 25 air
total.cumulative.consumption.PLOTTING[nrow(total.cumulative.consumption.PLOTTING) + 1,] <- 
  list("15", "25", "Heat stress", 15, 75, 6, "15°C water/25°C air", "Before")
##15 water and 30 air
total.cumulative.consumption.PLOTTING[nrow(total.cumulative.consumption.PLOTTING) + 1,] <- 
  list("15", "30", "Heat stress", 18, 15, 2, "15°C water/30°C air", "Before")
##20 water and 30 air
total.cumulative.consumption.PLOTTING[nrow(total.cumulative.consumption.PLOTTING) + 1,] <- 
  list("20", "30", "Heat stress", 13, 49, 7, "20°C water/30°C air", "Before")
##20 water and 25 air
total.cumulative.consumption.PLOTTING[nrow(total.cumulative.consumption.PLOTTING) + 1,] <- 
  list("20", "25", "Heat stress", 12, 88, 7, "20°C water/25°C air", "Before")
##20 water and 30 air
total.cumulative.consumption.PLOTTING[nrow(total.cumulative.consumption.PLOTTING) + 1,] <- 
  list("20", "20", "Heat stress", 10, 76, 14, "20°C water/20°C air", "Before")


#################################################
# cumulativeulative consumption per treatment per day
#################################################

#Making the plot
cumulativeconsumption_Plot_bytreatment <- total.cumulative.consumption.PLOTTING %>% 
  ggplot(aes(x = day.of.trial, y = total.cumulative.number.mussels, color = Treatment, linetype = changepoint)) +
  #Add the points and lines
  geom_line(size = 1) +
  scale_linetype_manual(name = "Change point",
                        values = c(1,4)) +
  geom_point(size = 1.5) +
  #Adjust the aesthetics 
  theme_classic() +
  scale_color_manual(name = "Treatment (°C)"
                     ,values = c("#9DC7C8", "#197BBD", "#05299E" , "#FC6DAB", "#FC814A", "#8D0801")) +
  scale_fill_manual(name = "Treatment (°C)"
                    ,values = c("#9DC7C8", "#197BBD", "#05299E" , "#FC6DAB", "#FC814A", "#8D0801")) +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8)) +
  guides(color = guide_legend(override.aes = list(size = 0.5))) +
  guides(fill = guide_legend(override.aes = list(size = 0.5))) +
  labs(x = "Day of experiment", y = "cumulativeulative number of mussels consumed per treatment") +
  theme(text = element_text(family = "A"),
        axis.title.x = element_text(vjust = -1, size = 12),
        axis.title.y = element_text(vjust = 2, size = 12),
        strip.text = element_text(size = 12),
        panel.background = element_rect(fill = NA, color = "black"),
        legend.box="vertical", 
        legend.margin=margin(t=0, unit ="cm")) +
  scale_x_continuous(breaks = seq(0, 24, by = 2)) +
  scale_y_continuous(breaks = seq(0, 200, by = 20))+
  #Make a rectangle for the recovery period
  annotate("rect", fill = "grey", alpha = 0.3, 
           xmin = 17, xmax = 24,
           ymin = -Inf, ymax = Inf) 


cumulativeconsumption_Plot_bytreatment



#################################################
# Total consumption per treatment per day
#################################################

#Making the plot
totalconsumption_Plot_bytreatment <- total.cumulative.consumption.PLOTTING %>% 
  ggplot(aes(x = day.of.trial, y = total.number.mussels, color = Treatment, linetype = changepoint)) +
  #Add the points and lines
  geom_line(size = 1) +
  scale_linetype_manual(name = "Change point",
                        values = c(1,4)) +
  geom_point(size = 1.5) +
  #Adjust the aesthetics 
  theme_classic() +
  scale_color_manual(name = "Treatment (°C)"
                     ,values = c("#9DC7C8", "#197BBD", "#05299E" , "#FC6DAB", "#FC814A", "#8D0801")) +
  scale_fill_manual(name = "Treatment (°C)"
                    ,values = c("#9DC7C8", "#197BBD", "#05299E" , "#FC6DAB", "#FC814A", "#8D0801")) +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8)) +
  guides(color = guide_legend(override.aes = list(size = 0.5))) +
  guides(fill = guide_legend(override.aes = list(size = 0.5))) +
  labs(x = "Day of experiment", y = "cumulativeulative number of mussels consumed per treatment") +
  theme(text = element_text(family = "A"),
        axis.title.x = element_text(vjust = -1, size = 12),
        axis.title.y = element_text(vjust = 2, size = 12),
        strip.text = element_text(size = 12),
        panel.background = element_rect(fill = NA, color = "black"),
        legend.box="vertical", 
        legend.margin=margin(t=0, unit ="cm")) +
  scale_x_continuous(breaks = seq(0, 24, by = 2)) +
  scale_y_continuous(breaks = seq(0, 200, by = 20))+
  #Make a rectangle for the recovery period
  annotate("rect", fill = "grey", alpha = 0.3, 
           xmin = 17, xmax = 24,
           ymin = -Inf, ymax = Inf) 


totalconsumption_Plot_bytreatment



##########################################################################
# Difference in cumulativeulative consumption between heat stress and recovery
## 15 water and 30 air (changepoint at day 18)
##########################################################################

## Looking at the last eight days of the heat stress period 
## and the first eight days of the recovery period 

##Set up the data for plotting
heatrecovery_plotting <- master.number %>% 
  filter(trial.type != "Adjustment") %>% 
  filter(air.treatment == "30" & water.treatment == "15" ) %>% 
  mutate(trial.block = case_when(
    day.of.trial %in% 1:8 ~ "Day 1-8",
    day.of.trial %in% 9:16 ~ "Day 9-16",
    day.of.trial %in% 17:24 ~ "Day 17-24",
    TRUE ~ "other"
  )) %>% 
  group_by(star.ID, trial.block, trial.type, water.treatment, air.treatment) %>% 
  summarise(total.number.mussels = sum(mussels.consumed)) 


#######################################
# Test for a difference in means
#######################################

#Start with an ANOVA
aov.changepoint <- aov(total.number.mussels ~ trial.block,
                       data = heatrecovery_plotting)

summary(aov.changepoint)

#Check homogeneity of variance assumption
leveneTest(total.number.mussels ~ trial.block,
           data = heatrecovery_plotting) #Variance is EQUAL

#Check normality assumption
##Extract residuals first
aov.changepoint_residuals <- residuals(object = aov.changepoint)
#Run shapiro wilk
shapiro.test(x = aov.changepoint_residuals) #NOT NORMAL 

##Kruskal-Wallis rank sum test (non-parametric alternative)

kruskal.test(total.number.mussels ~ trial.block,
             data = heatrecovery_plotting) # p = 0.01
#Post-hoc test
library(FSA)
dunnTest(total.number.mussels ~ trial.block,
         data = heatrecovery_plotting)

# Day 1-8 versus Day 17-24 (*)
# Day 1-8 versus Day 9-16 
# Day 17-24 versus Day 9-16 (*)

##########
# PLOT
##########

#Set up the data
heatrecovery_plotting$trial.block <- factor(heatrecovery_plotting$trial.block,
                                            ordered = TRUE,
                                            levels = c("Day 1-8", "Day 9-16", "Day 17-24"))

heatrecovery_plotting <- heatrecovery_plotting %>% 
mutate(Treatment = case_when(
  water.treatment =="15" & air.treatment=="20" ~ "15°C water/20°C air",
  water.treatment =="15" & air.treatment=="25" ~ "15°C water/25°C air", 
  water.treatment =="15" & air.treatment=="30" ~ "15°C water/30°C air", 
  water.treatment =="20" & air.treatment=="20" ~ "20°C water/20°C air",
  water.treatment =="20" & air.treatment=="25" ~ "20°C water/25°C air", 
  water.treatment =="20" & air.treatment=="30" ~ "20°C water/30°C air", 
  TRUE ~ "other"
))


## MAKE THE PLOT

changepoint.plot <- heatrecovery_plotting %>% 
  ggplot() +
  geom_boxplot(aes(x = trial.block, y = total.number.mussels, 
                   fill = trial.block, color = Treatment),
               outlier.shape = NA) +
  scale_color_manual(name = "Treatment (°C)", 
                     values = c("#05299E")) +
  scale_fill_manual(name = "Trial block",
                    values = c("white", "white", "grey91"))+
  LW_theme +
  guides(fill = "none") +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.box="vertical", 
        legend.margin=margin(t=0, unit ="cm")) +
  labs(x = "Experimental Period", y = "Total number of mussels consumed per individual") +
  scale_y_continuous(breaks = seq(0, 5, by = 1)) +
 # geom_text(label = "a", x = 1, y = 1.5, size = 5, family = "A") +
#  geom_text(label = "a", x = 2, y = 1.5, size = 5, family = "A") +
  geom_text(label = "*", x = 3, y = 4.5, size = 5, family = "A") 

changepoint.plot

  
#######################
# MS FIGURE
#######################

MS_cumulativefeeding_exp1 <- cumulativeconsumption_Plot_bytreatment + changepoint.plot +
  plot_layout(ncol = 2, widths = c(1, 0.5)) +
  plot_annotation(tag_levels = "a") 


MS_cumulativefeeding_exp1

#png("MS Figures/Fig.4_CumulativeMusselsConsumed_Exp1.png", width = 9, height = 7, units = "in", res = 600)
#MS_cumulativefeeding_exp1
#dev.off()
  

########################
# DATA EXPLANATION
########################

##CHANGEPOING ANALYSIS

# changepoint analysis was conducted on both cumulative consumption and total consumption data
## Both methods found changepoints on the same days
## Cumulative consumption was chosen as a visual because its more readable

# A Kruskal Wallace and Dunn's post-hoc was used to detect differences between 8-day blocks of the trial period
## A sign. diff (p<0.5) was detected between day1-8 versus day17-24
## A sign. diff (p<0.5) was detected between day9-16 versus day17-24
  
  
  

