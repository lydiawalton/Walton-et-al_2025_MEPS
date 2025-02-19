######################################################################################i
# Title: Investigating the impact of contrasting air and seawater temperatures on the 
#        mortality of juvenile Pisaster ochraceus 
#
# Author: Lydia Walton
# Last updated: 13-Dec-2024
######################################################################################i

#dev.off() #Clear open plots
rm(list = ls(all=T)) #Clear environment

#Set working directory
library(here)
#Packages for data manipulation
library(dplyr)
library(tidyr)
library(tibble)
#Packages for survival analysis
library(survival)
library(ranger)
library(survminer)
#Packages for data visualization
library(ggplot2)
library(ggfortify)
library(patchwork)
library(extrafont)

#----------------------------------Load data ----

##Time - survival time in days
##Status - censoring status; 1 = censored, 2 = dead

PO_mortality_June2023 <- read.csv("2023 Data/Pisaster_Mortality-data_June2023.csv")

PO_mortality_2023vs2024 <- read.csv("2023 Data/Pisaster_Mortality-data_2023versus2024.csv")

#---------------------------------Create theme for plotting ----

#Lydia's theme
LW_theme <- theme_classic() +
  theme(axis.title.x = element_text(vjust = -1, size = 12),
    axis.title.y = element_text(vjust = 2, size = 11),
    strip.text = element_text(size = 12),
    legend.text = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    panel.background = element_rect(fill = NA, color = "black"))

#---------------------------------Get animal silhouettes from phylopic ----
library(rphylopic)

seastar <- pick_phylopic(name = "Asteriidae", n = 1)

# PhyloPic has a lot of contributors and we should acknowledge 
# their work. You can get data about images using get_attribution

# Get valid uuid
###star.cont <- get_uuid(name = "Asteriidae", n=1)

# Get attribution data for uuid
####get_attribution(uuid = star.cont)


#######################################i
# Kaplan-Meier Survival Analysis ----
#######################################i

##--------------------------------------JUNE 2023 ----

#Compute the survival curve (Kaplan-Meier estimate)
km.fit.June2023 <- survfit(Surv(time, status) ~ water.treatment + air.treatment, 
                           data = PO_mortality_June2023)
print(km.fit.June2023)

# Summary of survival curves
summary(km.fit.June2023)
# Access to the sort summary table
summary(km.fit.June2023)$table

#Access the values returned by survfit()
mortality.PO.June2023 <- data.frame(time = km.fit.June2023$time,
                          n.risk = km.fit.June2023$n.risk,
                          n.event = km.fit.June2023$n.event,
                          n.censor = km.fit.June2023$n.censor,
                          surv = km.fit.June2023$surv,
                          upper = km.fit.June2023$upper,
                          lower = km.fit.June2023$lower)
head(mortality.PO.June2023)

#---------------------------Plot the survival curve (Kaplan Meier curve)
PO.mortality.June2023.Figure <- ggsurvplot(km.fit.June2023, data = PO_mortality_June2023,
                          # conf.int = TRUE, 
                          #pval = TRUE,
                          legend.labs = c("15°C water/20°C air", "15°C water/25°C air",
                                          "15°C water/30°C air", "20°C water/20°C air",
                                          "20°C water/25°C air", "20°C water/30°C air"),
                          legend = c(0.4,0.25),
                          palette = c("#9DC7C8", "#197BBD", "#05299E" , "#FC6DAB", "#FC814A", "#8D0801"),
                          legend.title = "Treatment",
                          break.time.by = 1,
                          xlab = "Time (days)")
# risk.table = TRUE,
# tables.theme = theme_cleantable()

##You can edit the KM curve in ggplot but you have to specify $plot
PO.mortality.June2023.Figure$plot <- PO.mortality.June2023.Figure$plot +
  #Add theme to plot
  ggplot2::theme(#text = element_text(family = "A"),
    axis.title.x = element_text(vjust = -1, size = 12),
    axis.title.y = element_text(vjust = 3, size = 12),
    axis.text = element_text(size = 8),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12)) +
  ggplot2::ylim(0.5,1.00) +
 # guides(color = guide_legend(nrow = 2)) +
  #Add p-value
  geom_text(label = "p-value < 0.01", x = 0.75, y = 0.52,
            size = 4) +
  #Add animal outline
  add_phylopic(x = 0.75, y = 0.58, img = seastar, alpha = 1, ysize = 0.1)

PO.mortality.June2023.Figure #Look at the plot

##### MS Figure ----

#png("MS Figures/Fig2_Mortality_June2023.png", width = 10, height = 9, units = "in", res = 600)
#PO.mortality.June2023.Figure
#dev.off()

##Cox proportional hazards model 

cox.PO.June2023 <- coxph(Surv(time, status)~water.treatment*air.treatment, 
                data = PO_mortality_June2023)

summary(cox.PO.June2023) #Significant differences between treatments

##Just looking at the 30air treatments
PO.30air.June2023 <- PO_mortality_June2023 %>% 
  filter(air.treatment == "30")

PO.cox.30air <- coxph(Surv(time, status)~water.treatment*air.treatment, 
                               data = PO.30air.June2023)

summary(PO.cox.30air) #No difference between 30air treatments


##Cox model w/out 30air and 15water

PO.cox.no30air15water <- PO_mortality_June2023 %>% 
  filter(water.treatment != "15" | air.treatment != "30")

cox.PO.no30air15water <- coxph(Surv(time, status)~water.treatment*air.treatment, 
                               data = PO.cox.no30air15water)

summary(cox.PO.no30air15water) # p < 0.05 (likelihood ratio) and p < 0.05 (logrank test)
###30air and 20water is different from other treatments

#Cox model w/out 30air and 20water

PO.cox.no30air20water <- PO_mortality_June2023 %>% 
  filter(water.treatment != "20" | air.treatment != "30")

cox.PO.no30air20water <- coxph(Surv(time, status)~water.treatment*air.treatment, 
                               data = PO.cox.no30air20water)

summary(cox.PO.no30air20water) # p = 0.005 (likelihood ratio) and p = 0.01 (logrank test)
###30 air and 15 water is different from other treatments


##--------------------------------------2023 versus 2024 ----

#Set the months to a specific order
PO_mortality_2023vs2024$experiment <- factor(PO_mortality_2023vs2024$experiment, order = TRUE,
                                             levels = c("June 2023", "August 2023", "May 2024"))

#Compute the survival curve (Kaplan-Meier estimate)
km.fit.2023vs2024 <- survfit(Surv(time, status) ~ experiment, 
                           data = PO_mortality_2023vs2024)
print(km.fit.2023vs2024)

# Summary of survival curves
summary(km.fit.2023vs2024)
# Access to the sort summary table
summary(km.fit.2023vs2024)$table

#Access the values returned by survfit()
mortality.PO.2023vs2024 <- data.frame(time = km.fit.2023vs2024$time,
                                    n.risk = km.fit.2023vs2024$n.risk,
                                    n.event = km.fit.2023vs2024$n.event,
                                    n.censor = km.fit.2023vs2024$n.censor,
                                    surv = km.fit.2023vs2024$surv,
                                    upper = km.fit.2023vs2024$upper,
                                    lower = km.fit.2023vs2024$lower)
head(mortality.PO.2023vs2024)

#---------------------------Plot the survival curve (Kaplan Meier curve)

#Make the baseplot using ggsurvplot
PO.mortality.2023vs2024.Figure <- ggsurvplot(km.fit.2023vs2024, data = PO_mortality_2023vs2024,
                                           # conf.int = TRUE, 
                                           #pval = TRUE,
                                           legend.labs = c("June 2023", "August 2023", "May 2024"),
                                           legend = c(0.4,0.25),
                                           palette = c("#26547C", "#53131E", "#E56399"),
                                           legend.title = "Experiment",
                                           break.time.by = 1,
                                           xlab = "Time (days)")
# risk.table = TRUE,
# tables.theme = theme_cleantable()

##You can edit the KM curve in ggplot but you have to specify $plot
PO.mortality.2023vs2024.Figure$plot <- PO.mortality.2023vs2024.Figure$plot +
  #Add theme to plot
  ggplot2::theme(text = element_text(family = "A"),
                 axis.title.x = element_text(vjust = -1, size = 12),
                 axis.title.y = element_text(vjust = 3, size = 12),
                 axis.text = element_text(size = 8),
                 panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 axis.line = element_line(colour = "black"),
                 legend.position = "bottom") +
 # guides(color = guide_legend(nrow = 2)) +
  #Add p-value
  geom_text(label = "p-value < 0.001", x = 1, y = 0.1,
            size = 4, family = "A") +
  #Add animal outline
  add_phylopic(x = 0.75, y = 0.25, img = seastar, alpha = 1, ysize = 0.2)

PO.mortality.2023vs2024.Figure #Look at the plot

##---------------------------June 2023 (Experiment 1) vs May 2024 (Experiment 2) ----
PO_mortality_June2023vsMay2024 <- PO_mortality_2023vs2024 %>% 
  filter(experiment != "August 2023")

#Compute the survival curve (Kaplan-Meier estimate)
km.fit.June2023vsMay2024 <- survfit(Surv(time, status) ~ experiment, 
                             data = PO_mortality_June2023vsMay2024)
print(km.fit.June2023vsMay2024)

# Summary of survival curves
summary(km.fit.June2023vsMay2024)
# Access to the sort summary table
summary(km.fit.June2023vsMay2024)$table

#Access the values returned by survfit()
mortality.PO.June2023vsMay2024 <- data.frame(time = km.fit.June2023vsMay2024$time,
                                      n.risk = km.fit.June2023vsMay2024$n.risk,
                                      n.event = km.fit.June2023vsMay2024$n.event,
                                      n.censor = km.fit.June2023vsMay2024$n.censor,
                                      surv = km.fit.June2023vsMay2024$surv,
                                      upper = km.fit.June2023vsMay2024$upper,
                                      lower = km.fit.June2023vsMay2024$lower)
head(mortality.PO.June2023vsMay2024)


#Make the baseplot using ggsurvplot
PO.mortality.June2023vsMay2024.Figure <- ggsurvplot(km.fit.June2023vsMay2024, data = PO_mortality_June2023vsMay2024,
                                             # conf.int = TRUE, 
                                             #pval = TRUE,
                                             legend.labs = c("June 2023", "May 2024"),
                                             legend = c(0.4,0.25),
                                             palette = c("#60D394", "#DB7F8E"),
                                             legend.title = "Experiment",
                                             break.time.by = 1,
                                             xlab = "Time (days)")
# risk.table = TRUE,
# tables.theme = theme_cleantable()

##You can edit the KM curve in ggplot but you have to specify $plot
PO.mortality.June2023vsMay2024.Figure$plot <- PO.mortality.June2023vsMay2024.Figure$plot +
  #Add theme to plot
  ggplot2::theme(text = element_text(family = "A"),
                 axis.title.x = element_text(vjust = -1, size = 12),
                 axis.title.y = element_text(vjust = 3, size = 12),
                 axis.text = element_text(size = 8),
                 panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 axis.line = element_line(colour = "black"),
                 legend.position = "right",
                 legend.title = element_text(size = 12),
                 legend.text = element_text(size = 12)) +
  # guides(color = guide_legend(nrow = 2)) +
  #Add p-value
  geom_text(label = "p-value < 0.01", x = 0.75, y = 0.1,
            size = 4, family = "A") +
  #Add animal outline
  add_phylopic(x = 0.75, y = 0.25, img = seastar, alpha = 1, ysize = 0.2)

PO.mortality.June2023vsMay2024.Figure

##Cox proportional hazards model 

cox.PO.June2023vsMay2024 <- coxph(Surv(time, status)~experiment, 
                         data = PO_mortality_June2023vsMay2024)

summary(cox.PO.June2023vsMay2024) #Significant differences between experiments (p < 0.01)












