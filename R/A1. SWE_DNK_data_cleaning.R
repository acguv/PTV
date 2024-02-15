
# ---------------------------------------------------------------------------- #
# Title:   Population total variation
# Country: Sweden and Denmark by Education
# Author:  Gomez-Ugarte Ana C 
# ---------------------------------------------------------------------------- #

# Content:
#   0. Working directory, packages and functions
#   1. Import data
#   2. Clean data
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
library(HMDHFDplus)

# Working directory from safe environment
#setwd("(...)")

# Import functions
source("R/Functions_PTV.R")

# ---------------------------------------------------------------------------- #
#     1. Import data
# ---------------------------------------------------------------------------- #

# 1. Life tables by education with weights
# Data provided by the authors of Nemeth et al. (2021)
Edu_LifeTables_weights <- read.csv("Data/DNKSWE_Ltables_eduprop_1991-2015.csv")


# 2. Total population and deaths for Denmark and Sweden
# Data downloaded from the Human Mortality Database (HMD) using the r package HMDHFDplus
username = "gomezugartevalerio@demogr.mpg.de"
pwd = "Pako2611!"

Dnk_pop <- readHMDweb("DNK","Population5",username = username,password = pwd)
Swe_pop <- readHMDweb("SWE","Population5",username = username,password = pwd) 


# ---------------------------------------------------------------------------- #
#     2. Clean data
# ---------------------------------------------------------------------------- #

# Clean and the arrange the population data 
Dnk_pop <- Dnk_pop |> 
  filter(Year %in% c(1991:2015)) |> 
  mutate(Country = "Denmark", 
         Period = as.numeric(cut_number(Year,5))) |> 
  mutate(Age = case_when(Age >= 90 ~ 90,
                         TRUE ~ as.numeric(Age))) |> 
  select(Country, Year, Age, Period, Female1, Male1, Total1) |>
  group_by(Country, Year, Age, Period) |> 
  summarise(Female1 = sum(Female1),
            Male1 = sum(Male1),
            Total1 = sum(Total1)) |> 
  group_by(Country, Period, Age) |> 
  summarise(Female1 = sum(Female1),
            Male1 = sum(Male1),
            Total1 = sum(Total1)) |>
  filter(Age >= 30) |>
  mutate(Period = case_when(Period == 1 ~ "1991-1995",
                            Period == 2 ~ "1996-2000",
                            Period == 3 ~ "2001-2005",
                            Period == 4 ~ "2006-2010",
                            Period == 5 ~ "2011-2015")) |>
  select(-Total1) |>
  rename(Female = Female1, Male = Male1) |>
  gather(Sex, population, Female:Male)


Swe_pop <- Swe_pop |> 
  filter(Year %in% c(1991:2015)) |> 
  mutate(Country = "Sweden", 
         Period = as.numeric(cut_number(Year,5))) |> 
  mutate(Age = case_when(Age >= 90 ~ 90,
                         TRUE ~ as.numeric(Age))) |> 
  select(Country, Year, Age, Period, Female1, Male1, Total1) |>
  group_by(Country, Year, Age, Period) |> 
  summarise(Female1 = sum(Female1),
            Male1 = sum(Male1),
            Total1 = sum(Total1)) |> 
  group_by(Country, Period, Age) |> 
  summarise(Female1 = sum(Female1),
            Male1 = sum(Male1),
            Total1 = sum(Total1)) |>
  filter(Age >= 30) |>
  mutate(Period = case_when(Period == 1 ~ "1991-1995",
                            Period == 2 ~ "1996-2000",
                            Period == 3 ~ "2001-2005",
                            Period == 4 ~ "2006-2010",
                            Period == 5 ~ "2011-2015")) |>
  select(-Total1) |>
  rename(Female = Female1, Male = Male1) |>
  gather(Sex, population, Female:Male)

# Merge data from Denmark and Sweden on population
Swe_Den_pop <- rbind(Dnk_pop, Swe_pop)


# Create life table of the total population by aggregating the mortality rates of the subgroups
mx_tot <- Edu_LifeTables_weights %>%
  select(Country, Period, Sex, Education, EduProp, Age, mx, ax) %>%
  left_join(Swe_Den_pop, by = c("Country", "Period", "Sex", "Age")) %>%
  mutate(EduPop = EduProp*population) %>%
  group_by(Country, Period, Sex, Age) %>%
  summarise(mx = sum(mx*EduPop)/sum(EduPop)) %>%
  arrange(Country, Period, Sex, Age) %>%
  group_by(Country, Period, Sex) %>%  
  mutate(dx = lifetable.mx(seq(30,90,5), mx, ax=NULL)$dx,
         ax = lifetable.mx(seq(30,90,5), mx, ax=NULL)$ax,
         ex = lifetable.mx(seq(30,90,5), mx, ax=NULL)$ex,
         mx = lifetable.mx(seq(30,90,5), mx, ax=NULL)$mx,
         lx = lifetable.mx(seq(30,90,5), mx, ax=NULL)$lx) %>%
  mutate(Education = "Total") 

# Merge total population life table to the one of all the groups
lt_edu <- Edu_LifeTables_weights %>%
  left_join(Swe_Den_pop, by = c("Country", "Period", "Sex", "Age")) %>%
  mutate(EduPop = EduProp*population) %>%
  select(-population) %>%
  bind_rows(mx_tot)  %>%
  mutate(AgeGrp = case_when(Age == 30 ~ "30-34",
                            Age == 35 ~ "35-39",
                            Age == 40 ~ "40-44",
                            Age == 45 ~ "45-49",
                            Age == 50 ~ "50-54",
                            Age == 55 ~ "55-59",
                            Age == 60 ~ "60-64",
                            Age == 65 ~ "65-69",
                            Age == 70 ~ "70-74",
                            Age == 75 ~ "75-79",
                            Age == 80 ~ "80-84",
                            Age == 85 ~ "85-89",
                            Age == 90 ~ "90+",
                            TRUE ~ "")) %>%
  group_by(Country, Sex, Education, Period) %>% 
  mutate(dx.d = dx/sum(dx)) %>%
  ungroup() %>%
  mutate(Education = fct_relevel(Education, "Low", "Middle", "High", "Total")) %>%
  arrange(Country, Period, Sex, Age, Education, Period) 

# Graph of mortality rates by Education level
lt_edu %>%
  filter(Sex == "Male" & Country == "Denmark" & Period == "2011-2015") %>%
  ggplot(aes(x=Age, y = log(mx), group = Education, color = Education))+
  geom_line()

# Estimate population weights and population relative ranks of the education groups
Edu_ranks <- Edu_LifeTables_weights %>%
  select(Country, Period, Sex, Education, EduProp, Age) %>%
  left_join(Swe_Den_pop, by = c("Country", "Period", "Sex", "Age")) %>%
  mutate(EduPop = EduProp*population) %>%
  group_by(Country, Period, Sex, Education) %>%
  summarise(TotEduPop = sum(EduPop)) %>%
  group_by(Country, Period, Sex) %>%
  mutate(Edu_weights = TotEduPop/sum(TotEduPop)) %>%
  select(-TotEduPop) %>%
  spread(Education, Edu_weights) %>%
  mutate(Low1 = 0.5*Low + Middle + High,
         Middle1 = 0.5*Middle + High,
         High1 = 0.5*High) |>
  select(Country, Period, Sex, Low1, Middle1, High1) %>%
  rename(Low = Low1, Middle = Middle1, High = High1) %>%
  gather(Education, Edu_ranks, Low:High) %>%
  mutate(Education = fct_relevel(Education, "Low", "Middle", "High")) %>%
  arrange(Country, Period, Sex, Education)

Edu_weights <- Edu_LifeTables_weights |>
  select(Country, Period, Sex, Education, EduProp, Age) |>
  left_join(Swe_Den_pop, by = c("Country", "Period", "Sex", "Age")) |>
  mutate(EduPop = EduProp*population) |>
  group_by(Country, Period, Sex, Education) |>
  summarise(TotEduPop = sum(EduPop)) |>
  group_by(Country, Period, Sex) |>
  mutate(Edu_weights = TotEduPop/sum(TotEduPop)) |>
  mutate(Education = fct_relevel(Education, "Low", "Middle", "High")) %>%
  select(-TotEduPop) %>%
  arrange(Country, Period, Sex, Education)

# Graph of population weights by country and sex
Edu_weights %>%
  mutate(Education = fct_relevel(Education, "Low", "Middle", "High")) %>%
  ggplot(aes(x = Period, y = Edu_weights, fill = Education)) + 
  geom_bar(position = "stack", stat = "identity") +
  facet_grid(Sex ~ Country)


# ---------------------------------------------------------------------------- #
#     3. Save results
# ---------------------------------------------------------------------------- #

# save(lt_edu, file = "Data/SWE_DNK_LT_tot.RData")
# save(Edu_ranks, file = "Data/SWE_DNK_Edu_ranks.RData")
# save(Edu_weights, file = "Data/SWE_DNK_Edu_weigths.RData")

