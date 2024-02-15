
# ---------------------------------------------------------------------------- #
# Title:   Population total variation
# Country: Denmark and Sweden by education level
# Author:  Gomez-Ugarte Ana C 
# ---------------------------------------------------------------------------- #

# Content:
#   0. Working directory, package and functions
#   1. Read the data
#   2. Main figures

# ---------------------------------------------------------------------------- #
#     0. Working directory, package and functions
# ---------------------------------------------------------------------------- #

# To clear everything in R, before start the analysis
rm(list = ls()) 

# Libraries
library(ggplot2)
library(dplyr)
library(broom)
library(ggpubr)
library(tidyverse)
library(scales)

# Working directory from safe environment
#setwd("(...)")

# Import functions
source("R/Functions_PTV.R")

# ---------------------------------------------------------------------------- #
#     1. Read and prepare data
# ---------------------------------------------------------------------------- #

# Read data
load("Data/SWE_DNK_LT_tot.RData")
load("Data/SWE_DNK_Edu_ranks.RData")
load("Data/SWE_DNK_Edu_weigths.RData")

who_std <- read.csv("Data/WHO_std_5_age.csv")

# ---------------------------------------------------------------------------- #
#     2. Main figures
# ---------------------------------------------------------------------------- #

# ------Range and ratio of life expectancy------# 
ex_range_edu <- lt_edu %>%
  filter(AgeGrp == "30-34") %>%
  select(Country, Sex, Education, Period, ex) %>%
  spread(Education, ex) %>%
  # Range (absolute measure)
  mutate(gap = High-Low) %>%
  # Ratio (relative measure)
  mutate(ratio = High/Low) %>%
  select(Country, Sex, Period, gap, ratio)

# ------Range and ratio of lifespan variation------# 

# Vector with ages, to be used inside sdv_func 
Age <- seq(30,90,5)

sdv_range_edu <- lt_edu %>%
  arrange(Country, Sex, Education, Period, AgeGrp) %>%
  group_by(Country, Sex, Education, Period) %>%
  group_modify(~ sdv_func(.)) %>%
  ungroup()  %>%
  select(Country, Sex, Education, Period, `S[1]`) %>%
  unique() %>%
  spread(Education, `S[1]`) %>%
  # Range (absolute measure)
  mutate(gap = Low-High) %>%
  # Ratio (relative measure)
  mutate(ratio = Low/High) %>%
  select(Country, Sex, Period, gap, ratio) 

# ------Slope/Relative index of inequality for life expectancy------# 

who_std <- who_std[c(8:20),] %>%
  mutate(Prop = Prop/sum(Prop))
edu_lt <- merge(lt_edu, who_std, by = "AgeGrp")

# With real population shares by education level
sii_edu <- edu_lt %>%
  mutate(mx_w = mx*Prop) %>%
  group_by(Country, Sex, Education, Period) %>% 
  mutate(asmr = sum(mx_w)*100000)  %>%
  select(Country, Sex, Education, Period, asmr) %>%
  unique() %>%
  mutate(cat_edu = case_when(Education == "Low" ~ 1, 
                             Education == "Middle" ~ 2,
                             Education == "High" ~ 3)) %>%
  group_by(Country, Sex, Period) %>%
  arrange(Country, Sex, Period, cat_edu) %>%
  left_join(Edu_ranks, by = c("Country", "Period", "Sex", "Education")) %>%
  do(tidy(lm(asmr ~ Edu_ranks, .))) %>%
  select(-std.error, -statistic, -p.value) %>%
  spread(term, estimate) %>%
  # Slope index of inequality
  mutate(rii = (`(Intercept)` + Edu_ranks)/(`(Intercept)`)) %>%
  # Relative index of inequality
  mutate(sii = Edu_ranks) %>%
  select(Country, Sex, Period, rii, sii) %>%
  gather(measure, value, rii:sii) 

# ------Population Total Variation (PTV)------# 
ptv_edu <- lt_edu %>%
  arrange(Country, Sex, Period, Education) %>%
  left_join(Edu_weights, by = c("Country", "Period", "Sex", "Education")) %>%
  group_by(Country, Period, Sex, Age) %>% 
  mutate(M = dx.d[Education == "Total"]) %>%
  mutate(diff = pmax(dx.d, M)-pmin(dx.d, M)) %>%
  group_by(Country, Period, Sex, Education) %>%  
  mutate(ptv = 1/2*sum(diff)) %>%
  select(Country, Period, Sex, Education, Edu_weights, ptv) %>%
  unique() %>%
  group_by(Country, Period, Sex) %>%  
  filter(Education != "Total") %>%
  mutate(NPTV = sum(Edu_weights*ptv)) %>%
  select(Period, Sex, NPTV) %>%
  unique() 

# ---------------------------------------------------------------------------- #
#     3. Save results
# ---------------------------------------------------------------------------- #

# save(ex_range_edu, sdv_range_edu, sii_edu, ptv_edu,
#      file = "Results/SWE_DNK_indicators.RData")
