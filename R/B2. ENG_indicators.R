
# ---------------------------------------------------------------------------- #
# Title:   Population total variation
# Country: England by deprivation deciles
# ---------------------------------------------------------------------------- #

# Content:
#   0. Working directory, packages and functions
#   1. Read the data
#   2. Main figures
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

# Working directory from safe environment
#setwd("(...)")

# Import functions
source("R/Functions_PTV.R")

# ---------------------------------------------------------------------------- #
#     1. Read and prepare data
# ---------------------------------------------------------------------------- #

# 1. Life tables by decile Read data
load("Data/ENG_LT_TOT.RData")

# 2. Standard population
who_std <- read.csv("Data/WHO_std_single_age_90+.csv")
names(who_std)[1] <- "Age"

# ---------------------------------------------------------------------------- #
#     2. Main figures
# ---------------------------------------------------------------------------- #

# ------Range and ratio of life expectancy------# 
ex_range_dpi <- lt_dpi %>%
  filter(Age == 0 & Decile != "Total") %>%
  select(Period, Sex, Decile, ex) %>%
  spread(Decile, ex) %>%
  # Range (absolute measure)
  mutate(gap = `Decile 10`-`Decile 1`) %>%
  # Ration (relative measure)
  mutate(ratio = `Decile 10`/`Decile 1`) %>%
  select(Period,Sex, gap, ratio) %>%
  arrange(Sex, Period)

# ------Range and ratio of lifespan variation------# 

# Vector with ages, to be used inside sdv_func 
Age <- c(0:90)

sdv_range_dpi <- lt_dpi %>%
  filter(Decile != "Total") %>%
  group_by(Period, Sex, Decile) %>%
  group_modify(~ sdv_func(.)) %>%
  ungroup()  %>%
  select(Sex, Decile, Period, `S[1]`) %>%
  unique() %>%
  spread(Decile, `S[1]`) %>%
  # Range (absolute measure)
  mutate(gap = `Decile 1`-`Decile 10`) %>%
  # Ratio (relative measure)
  mutate(ratio = `Decile 1`/`Decile 10`) %>%
  select(Period,Sex, gap, ratio) %>%
  arrange(Sex, Period)

# ------Slope/Relative index of inequality for life expectancy------# 

# Merge life_tables by deprivation to the standard population
lt_dpi_std <- merge(lt_dpi, who_std, by = "Age")

# Population ranks, given that we are working with deciles, we can assume that we have 10% of the population in each decile. 
pop_rank <- seq(0.05, .95, .1)

sii_dpi <- lt_dpi_std %>%  
  mutate(mx_w = mx*Prop) %>%
  group_by(Period, Sex, Decile) %>% 
  mutate(asmr = sum(mx_w)*100000)  %>%
  select(Period, Sex, Decile, decile_num, asmr) %>%
  unique() %>%
  filter(Decile != "Total") %>%
  group_by(Sex, Period) %>%
  arrange(Sex, Period, -decile_num) %>%
  do(tidy(lm(asmr ~ pop_rank, .))) %>%
  select(-std.error, -statistic, -p.value) %>%
  spread(term, estimate) %>%
  mutate(rii = (`(Intercept)` + pop_rank)/(`(Intercept)`)) %>%
  mutate(sii = pop_rank) %>%
  select(Sex, Period, rii, sii) %>%
  gather(measure, value, rii:sii) 

# ------Population Total Variation (PTV)------# 
ptv_dpi <- lt_dpi %>%
  group_by(Year, Sex, Age) %>%  
  mutate(M = dx.d[Decile == "Total"]) %>%
  mutate(diff = pmax(dx.d, M)-pmin(dx.d, M)) %>%
  group_by(Period, Sex, Decile) %>%  
  mutate(ptv = 1/2*sum(diff)) %>%
  select(Period, Sex,Decile, ptv) %>%
  unique() %>%
  group_by(Period, Sex) %>%  
  mutate(NPTV = sum(1/10*ptv)) %>%
  select(Period, Sex, NPTV) %>%
  unique() %>%
  arrange(Sex, Period)


# ---------------------------------------------------------------------------- #
#     3. Save results
# ---------------------------------------------------------------------------- #
# save(ex_range_dpi, sdv_range_dpi, sii_dpi, ptv_dpi, file = "Results/ENG_indicators.RData")
