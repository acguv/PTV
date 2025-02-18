
# ---------------------------------------------------------------------------- #
# Paper:   Reassessing socioeconomic inequalities in mortality via distributional similarities
# Title:   Estimating measures of SES inequality in mortality
# Country: Denmark and Sweden by education level
# ---------------------------------------------------------------------------- #

# Content:
#   0. Working directory, package and functions
#   1. Read the data
#   2. Estimate measures
#   3. Save results

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
library(strat)

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
#     2. Estimate measures
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

# ------ Pairwise non-overlap index (Shi et al.)------#  
noi_pair_edu <- lt_edu %>%
  filter(Education != "Total") %>%
  left_join(Edu_weights, by = c("Country", "Period", "Sex", "Education")) %>%
  select(Period,  Country, Education, Sex, Age, dx.d, lx, Edu_weights, ex) %>%
  mutate(lx = lx/100000) %>%
  arrange(Period, Country, Sex, Education, Age) %>%
  group_by(Country, Period, Sex) %>%
  # filter(Period == "2006-08" & Sex == "Females") %>%
  mutate(noi_pair = noi_pair_func_dx(c(dx.d, Edu_weights), n = 3))

# ------Total non-overlap index (Shi et al)------# 
noi_edu <- lt_edu %>%
  arrange(Country, Sex, Period, Education) %>%
  group_by(Country, Period, Sex, Age) %>% 
  select(Country, Period, Sex,  Education,  Age, dx.d) %>%
  spread(Education, dx.d) %>%
  mutate(across("Low":"High" , function(x) apply(cbind(x, Total),1,max), .names = "max_{col}"),
         across("Low":"High" , function(x) apply(cbind(x, Total),1,min), .names = "min_{col}")) %>%
  select(Period, Sex, Age, starts_with("max"), starts_with("min")) %>%
  gather("serie", "value", "max_Low":"min_High") %>%
  mutate(limit = substr(serie, 1, 3),
         Education = substr(serie, 5, 13)) %>%
  group_by(Country, Period, Education, limit, Sex) %>%  
  summarise(add = sum(value)) %>%
  spread(limit, add) %>%
  group_by(Country, Period, Education, Sex) %>%  
  summarise(noi = 1-(sum(min)/sum(max))) %>%
  left_join(Edu_weights, by = c("Country", "Period", "Sex", "Education")) %>%
  group_by(Country, Period, Sex) %>% 
  mutate(w_all = sum(Edu_weights)) %>%
  mutate(cte_fact = 1-(Edu_weights/(1+w_all-Edu_weights))) %>%
  summarise(noi_t = sum(noi*Edu_weights)*(1/sum(Edu_weights*cte_fact))) %>%
  select(Country, Period, Sex, noi_t) %>%
  unique() %>%
  arrange(Country, Period, Sex)

# ------Stratification index (Zhou et al.)------# 
S_edu <- lt_edu %>%
  filter(Education != "Total") %>%
  arrange(Country, Sex, Period, Education) %>%
  group_by(Country, Period, Sex, Age) %>% 
  select(Country, Period, Sex,  Education,  Age, dx.d, EduPop) %>%
  group_by(Country, Period, Sex) %>%
  summarise(S = strat(Age, Education, weights = dx.d)$overall[1])

# ------Outsurvival probability (Vaupel et al.) for multiple populations------# 
OV_pair_edu <- lt_edu %>%
  filter(Education != "Total") %>%
  left_join(Edu_weights, by = c("Country", "Period", "Sex", "Education")) %>%
  select(Period,  Country, Education, Sex, Age, dx.d, lx, Edu_weights, ex) %>%
  mutate(lx = lx/100000) %>%
  arrange(Period, Country, Sex, Education, Age) %>%
  group_by(Country, Period, Sex) %>%
  # filter(Period == "2006-08" & Sex == "Females") %>%
  summarise(OV_pair = ov_pair_func_dx(dx.d, Edu_weights, lx, ex, n = 3))

# ------Population Total Variation (PTV) (for working paper)------# 
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

# save(ex_range_edu, sdv_range_edu, sii_edu, ptv_edu, noi_edu, noi_pair_edu, S_edu, OV_pair_edu,
#      file = "Results/SWE_DNK_indicators.RData")
