

# ---------------------------------------------------------------------------- #
# Title:   Population total variation
# Country: England
# ---------------------------------------------------------------------------- #

# Content:
#   0. Working directory, package and functions
#   1. Read the data
#   2. Decomposition
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
library(DemoDecomp)

# Working directory from safe environment
#setwd("(...)")

# Import functions
source("R/Functions_PTV.R")

# ---------------------------------------------------------------------------- #
#     1. Read and prepare data
# ---------------------------------------------------------------------------- #

# Load data
# 1. Life tables by decile Read data
load("Data/ENG_LT_TOT.RData")

# ---------------------------------------------------------------------------- #
#     2. Decomposition
# ---------------------------------------------------------------------------- #

#-------Decomposition of the PTV-----------#
 
# Function to get the the age-at-death-distribution in each group for a specific sex, country and period
get_mx = function(sex, year){
  mx <- lt_dpi %>%  
    filter(Decile != "Total") %>%
    select(Period, Year, Sex, Age, Decile, decile_num,mx) %>%
    arrange(Sex, Year, decile_num) %>% 
    filter(Sex == sex, Year == year) %>%
    .$mx
  return(mx)
}

# Function to get the populations shares in each group for a specific sex, country and period
get_weights = function(sex, year){
  wts <- lt_dpi %>%
    filter(Decile != "Total") %>%
    arrange(Sex, Year, Decile, decile_num) %>%
    filter(Sex == sex, Year == year) %>%
    .$Population
  return(wts)
}

# Getting the mortality rates and the population weights for the selected years
mx_f_07 <- get_mx(sex = "Females", year = "2007")
mx_f_15 <- get_mx(sex = "Females", year = "2015")

wts_f_07 <- get_weights(sex = "Females", year = "2007")
wts_f_15 <- get_weights(sex = "Females", year = "2015")

mx_w_f_07 <- c(mx_f_07, wts_f_07)
mx_w_f_15 <- c(mx_f_15, wts_f_15)

mx_m_07 <- get_mx(sex = "Males", year = "2007")
mx_m_15 <- get_mx(sex = "Males", year = "2015")

wts_m_07 <- get_weights(sex = "Males", year = "2007")
wts_m_15 <- get_weights(sex = "Males", year = "2015")

mx_w_m_07 <- c(mx_m_07, wts_m_07)
mx_w_m_15 <- c(mx_m_15, wts_m_15)

# Function to estimate the PTV given vectors with the age-at-death distributions (dx) and the population shares in each subgroup (w)
ptv_func <- function(mx_w, n = 10){
  w <- matrix(as.matrix(mx_w[(length(mx_w)/2+1):length(mx_w)]), ncol = n)
  w_dpi <- colSums(w)/sum(w)
  mx <- matrix(as.matrix(mx_w[1:(length(mx_w)/2)]), ncol = n)
  mx_tot <- rowSums(mx*w)/rowSums(w)
  dx_tot <- lifetable.mx(x = c(0:90), mx_tot, ax=NULL)$dx
  dx <- apply(mx, 2, function(mx) lifetable.mx(x=c(0:90), mx, ax=NULL)$dx)
  dx <- cbind(dx, dx_tot)
  dx.d <- dx/colSums(dx)
  diffs <- apply(dx.d,2,function(x) x-dx.d[,n+1])
  tv <- apply(diffs,2, function(x) sum(abs(x), na.rm = TRUE)/2)
  ptv <- sum(tv[-(n+1)]*w_dpi)
  ptv
}


# Decomposition by stepwise repalcement algorithm
decomp_f <- stepwise_replacement(ptv_func, mx_w_f_07, mx_w_f_15) 
decomp_m <- stepwise_replacement(ptv_func, mx_w_m_07, mx_w_m_15) 


# Function to reorder the results of the decomposition
contr_fun <- function(decomp_results, sex){
  contr <- data.frame(cbind(matrix(decomp_results[1:(length(decomp_results)/2)], ncol = 10), 
                            rowSums(matrix(decomp_results[((length(decomp_results)/2)+1):length(decomp_results)], ncol = 10))))
  names(contr) <- c(1:10, "Composition")
  rel_contr <- contr/sum(contr)
  names(rel_contr) <- c(1:10, "Composition")
  contr <- rbind(contr, rel_contr)
  contr %>%
    mutate(Sex = sex,
           Series = rep(c("Absolute", "Relative"), each = 91),
           Age = rep(c(0:90),2)) %>%
    gather("Decile", "Contribution", "1":"Composition") 
}

# Reordering the resutls of the decomposition
contr_uk_m <- contr_fun(decomp_m, "Males")
contr_uk_f <- contr_fun(decomp_f, "Females")

contr_uk_all_step <- rbind(contr_uk_f, contr_uk_m) 
  
# Graph of the decomposition by age and component
contr_uk_all_step %>%
  filter(Series == "Absolute") %>%
  mutate(Decile = fct_relevel(Decile, c("1", "2","3","4","5","6","7","8","9","10","Composition"))) %>%
  ggplot(aes(x= Age, y = Contribution, group = Decile, fill = Decile)) +
  geom_bar(stat = "identity") +
  facet_grid(~ Sex) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12), axis.title =  element_text(size = 13, face = "bold"), legend.text = element_text(size = 13), 
        plot.title = element_text(face = "bold"), legend.title = element_text(size = 13),
        axis.text.y = element_text(size = 12)) + xlab("Age") + ylab("Change in PTV") + 
  scale_fill_viridis_d() +
  theme(panel.grid.major.y = element_line(colour = "lightgrey")) +
  scale_x_continuous(breaks = seq(0, 90, by = 10))

ggsave("Graphs/ENG_decomp_ptv.pdf", width = 10.24, height =  4.64)

# Graph of the decomposition by component
contr_uk_all_step %>%
  filter(Series == "Absolute") %>%
  group_by(Sex, Decile) %>%
  summarise(Contribution = sum(Contribution)) %>%
  mutate(Decile = fct_relevel(Decile, c("1", "2","3","4","5","6","7","8","9","10","Composition"))) %>%
  ggplot(aes(x= Sex, y = Contribution, group = Decile, fill = Decile)) +
  geom_bar(stat = "identity", color = "black", position = "dodge") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12), axis.title =  element_text(size = 13, face = "bold"), legend.text = element_text(size = 13), 
        plot.title = element_text(face = "bold"), legend.title = element_text(size = 13),
        axis.text.y = element_text(size = 12)) + xlab("Sex") + ylab("Change in PTV") + 
  scale_fill_viridis_d() +
  theme(panel.grid.major.y = element_line(colour = "lightgrey")) 

#-------Decomposition of the SII-----------#

# Read standard population
who_std <- read.csv("Data/WHO_std_single_age_90+.csv")
names(who_std)[1] <- "Age"
who_std <- who_std[0:91, ]

# Function to estimate the PTV given vectors with the age-at-death distributions (dx) and the population shares in each subgroup (w)
sii_func <- function(mx_w, n = 10){
  w <- matrix(as.matrix(mx_w[(length(mx_w)/2+1):length(mx_w)]), ncol = n)
  w_dpi <- colSums(w)/sum(w)
  w_ranks <- seq(0.95,0.05,-0.1)
  mx <- matrix(as.matrix(mx_w[1:(length(mx_w)/2)]), ncol = n)
  asmr <- colSums(mx*who_std$Prop)*100000
  sii <- lm(asmr ~ w_ranks)$coefficients[2]
  sii
}

# Decomposition by stepwise replacement algorithm
decomp_sii_f <- stepwise_replacement(sii_func, mx_w_f_07, mx_w_f_15) 
decomp_sii_m <- stepwise_replacement(sii_func, mx_w_m_07, mx_w_m_15) 

# Reordering the results from the decomposition
contr_uk_m <- contr_fun(decomp_sii_m, "Males")
contr_uk_f <- contr_fun(decomp_sii_f, "Females")

contr_uk_all_sii_step <- rbind(contr_uk_f, contr_uk_m) 

# Graph of the decomposition by age and component
contr_uk_all_sii_step %>%
  filter(Series == "Absolute") %>%
  mutate(Decile = fct_relevel(Decile, c("1", "2","3","4","5","6","7","8","9","10","Composition"))) %>%
  ggplot(aes(x= Age, y = Contribution, group = Decile, fill = Decile)) +
  geom_bar(stat = "identity") +
  facet_grid(~ Sex) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12), axis.title =  element_text(size = 13, face = "bold"), legend.text = element_text(size = 13), 
        plot.title = element_text(face = "bold"), legend.title = element_text(size = 13),
        axis.text.y = element_text(size = 12)) + xlab("Sex") + ylab("Change in SII") + 
  scale_fill_viridis_d() +
  theme(panel.grid.major.y = element_line(colour = "lightgrey")) + 
  scale_x_continuous(breaks = seq(0, 90, by = 10))

ggsave("Graphs/ENG_decomp_sii.pdf", width = 10.24, height =  4.64)

# Graph of the decomposition by component
contr_uk_all_sii_step %>%
  filter(Series == "Absolute") %>%
  group_by(Sex, Decile) %>%
  summarise(Contribution = sum(Contribution)) %>%
  mutate(Decile = fct_relevel(Decile, c("1", "2","3","4","5","6","7","8","9","10","Composition"))) %>%
  ggplot(aes(x= Sex, y = Contribution, group = Decile, fill = Decile)) +
  geom_bar(stat = "identity", color = "black", position = "dodge") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12), axis.title =  element_text(size = 13, face = "bold"), legend.text = element_text(size = 13), 
        plot.title = element_text(face = "bold"), legend.title = element_text(size = 13),
        axis.text.y = element_text(size = 12)) + xlab("Sex") + ylab("Change in SII") + 
  scale_fill_viridis_d() +
  theme(panel.grid.major.y = element_line(colour = "lightgrey")) 

# ---------------------------------------------------------------------------- #
#     3. Save results
# ---------------------------------------------------------------------------- #
# save(contr_uk_all_step, contr_uk_all_sii_step, file = "Results/ENG_decomp_results_stepwise.RData")
