################################################################i
# Title: Moribundity surveys for Pisaster
#
# Author: Lydia Walton
# Last updated: 22-Jun-2025
################################################################i

# Setup ----
#dev.off
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
library(sjPlot)
#Statistics
library(glmmTMB)
library(DHARMa)

# Setup ----

#------------------------------------------Load data files

PO.moribundity <- read.csv("2023 data/Pisaster_MoribunditySurveys_June2023.csv")

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


## Statistical tests ----

### Chi-squared test ----

#### Sites ----

#Make contingency table
PO.moribundity.site <- PO.moribundity %>% 
  group_by(site) %>% 
  summarise(total.healthy = sum(total.healthy), total.moribund = sum(total.moribund))

library(tidyverse)
site.contable <- PO.moribundity.site %>% 
                    remove_rownames %>% 
                    column_to_rownames(var="site")

#Chi-square test
site.chisq <- chisq.test(site.contable)
site.chisq

#Extract observed counts
site.chisq$observed

#Extract expected counts
round(site.chisq$expected,2)

#Extract Pearson residuals
round(site.chisq$residuals,3)

#Visualize Pearson residuals
library(corrplot)
corrplot(site.chisq$residuals, is.cor = FALSE)

#Slight association between site and moribundity 
# - Eagle Bay has positive association with moribundity
# - Strawberry Point has positive association with healthy individuals

#We do have a low sample size at some sites
#Run a Monte-Carlo simulation
site.chisq.MC <- chisq.test(site.contable, simulate.p.value = TRUE, B = 10000)
site.chisq.MC


#### Month ----

#Make contingency table
PO.moribundity.month <- PO.moribundity %>% 
  group_by(month) %>% 
  summarise(total.healthy = sum(total.healthy), total.moribund = sum(total.moribund))

month.contable <- PO.moribundity.month %>% 
  remove_rownames %>% 
  column_to_rownames(var="month")

#Chi-square test
month.chisq <- chisq.test(month.contable)
month.chisq

#Extract observed counts
month.chisq$observed

#Extract expected counts
round(month.chisq$expected,2)

#Extract Pearson residuals
round(month.chisq$residuals,3)

#Visualize Pearson residuals
corrplot(month.chisq$residuals, is.cor = FALSE)

#Association between month and moribundity 
# - August has positive association with moribundity

#Month is not independent as sites were sampled multiple times so we violate that assumption of the chi-square test

#Run a Monte-Carlo simulation
month.chisq.MC <- chisq.test(month.contable, simulate.p.value = TRUE, B = 10000)
month.chisq.MC









