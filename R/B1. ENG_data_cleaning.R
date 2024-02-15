
# ---------------------------------------------------------------------------- #
# Title:   Population total variation
# Country: England
# Author:  Gomez-Ana C 
# ---------------------------------------------------------------------------- #

# Content:
#   1. Load packages and import data
#   2. Clean data
#   3. Save results

rm(list = ls())

# ---------------------------------------------------------------------------- #
#     1. Load packages and import data
# ---------------------------------------------------------------------------- #
library(ggplot2)
library(dplyr)
library(broom)
library(ggpubr)
library(scales)
library(ungroup)
library(tidyr)
library(zoo)


# 1. Population in England by sex and age
# pop <- read.csv("ENG/Population_Eng.csv")
pop <- read.csv("Data/Eng_Population_Deciles.csv")

# 3. UK_deprivation_deciles.csv: Life tables by deciles 
# Download data from the supplementary information from the Office of National Statistics 
# Available at: https://www.ons.gov.uk/peoplepopulationandcommunity/healthandsocialcare/healthandlifeexpectancies/adhocs/009149lifetablebysingleyearofagesexanddeprivationdecilesinenglandbetween2006to2008and2014to2016 
lt_eng_dep <- read.csv("Data/ENG_LT_deprivation_deciles.csv", header = TRUE, skip = 6)

# Import functions
source("R/Functions_PTV.R")

# ---------------------------------------------------------------------------- #
#     1. Clean data
# ---------------------------------------------------------------------------- #

# Clean data from life tables by decile and merge with the age-at-death distribution from the total population
lt_imd <- lt_eng_dep[-c(16381:16400),]

names(lt_imd) <- c("Period", "Decile", "Sex",  "Age", "mx", "ax", "qx", "lx", "dx", "Lx", "Tx", "ex", "ex_low", "ex_high")

lt_imd <- lt_imd %>%
  mutate(age = case_when(Age == "<1" ~ 0,
                         Age == "90+" ~ 90,
                         TRUE ~ as.numeric(Age)))%>%
  mutate(Decile = fct_relevel(Decile, "Decile 1","Decile 2","Decile 3","Decile 4","Decile 5","Decile 6","Decile 7","Decile 8",
                              "Decile 9","Decile 10")) %>%
  mutate(decile_num = as.numeric(substr(Decile, 8, 9))) %>%
  select(Period, Sex, age, Decile, decile_num, dx, ax, ex, lx, mx) %>%
  rename(Age = age) %>%
  mutate(Year = case_when(Period == "2006-08" ~ 2007,
                          Period == "2007-09" ~ 2008,
                          Period == "2008-10" ~ 2009,
                          Period == "2009-11" ~ 2010,
                          Period == "2010-12" ~ 2011,
                          Period == "2011-13" ~ 2012,
                          Period == "2012-14" ~ 2013,
                          Period == "2013-15" ~ 2014,
                          Period == "2014-16" ~ 2015,
                          TRUE ~ NA_real_)) %>%
  arrange(Year, Sex, decile_num, Age) 


mx_tot <- lt_imd %>%
  left_join(pop, by = c("Year", "Sex", "decile_num", "Age")) %>%
  group_by(Year, Sex, Age) %>%
  summarise(mx = sum(mx*Population)/sum(Population)) %>%
  arrange(Year, Sex, Age) %>%
  group_by(Year, Sex) %>%  
  mutate(dx = lifetable.mx(x = Age, mx = mx, ax=NULL)$dx,
         ax = lifetable.mx(x = Age, mx = mx, ax=NULL)$ax,
         ex = lifetable.mx(x = Age, mx = mx, ax=NULL)$ex,
         mx = lifetable.mx(x = Age, mx = mx, ax=NULL)$mx,
         lx = lifetable.mx(x = Age, mx = mx, ax=NULL)$lx) %>%
  mutate(Decile = "Total") 
  
lt_dpi <- lt_imd %>%
  left_join(pop, by = c("Year", "Sex", "decile_num", "Age")) %>%
  arrange(Year, Sex, Decile, Age) %>% 
  bind_rows(mx_tot) %>%
  group_by(Year, Sex, Decile) %>% 
  mutate(dx.d = dx/sum(dx)) %>%
  ungroup() %>%
  mutate(decile_num = case_when(is.na(decile_num) ~ 0,
                                TRUE ~ decile_num)) %>%
  filter(Year > 2006 & Year < 2016) %>%
  mutate(Period = case_when(Year == 2007 ~ "2006-08",
                            Year == 2008 ~ "2007-09",
                            Year == 2009 ~ "2008-10",
                            Year == 2010 ~ "2009-11",
                            Year == 2011 ~ "2010-12",
                            Year == 2012 ~ "2011-13",
                            Year == 2013 ~ "2012-14",
                            Year == 2014 ~ "2013-15",
                            Year == 2015 ~ "2014-16",
                          TRUE ~ "")) %>%
  arrange(Year, Sex, decile_num, Age) %>%
  mutate(Decile = fct_relevel(Decile, "Decile 1","Decile 2","Decile 3","Decile 4","Decile 5","Decile 6","Decile 7","Decile 8",
                              "Decile 9","Decile 10", "Total")) 
  
# ---------------------------------------------------------------------------- #
#     3. Save results
# ---------------------------------------------------------------------------- #

# save(lt_dpi, file = "Data/ENG_LT_TOT.RData")

