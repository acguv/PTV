
# ---------------------------------------------------------------------------- #
# Title:   Population total variation
# Country: England by deprivation deciles
# Author:  Gomez-Ugarte Ana C 
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
load("Results/ENG_indicators.RData")

# 2. Standard population
who_std <- read.csv("Data/WHO_std_single_age_90+.csv")
names(who_std)[1] <- "Age"

# ---------------------------------------------------------------------------- #
#     2. Main figures
# ---------------------------------------------------------------------------- #

# ------Range and ratio of life expectancy------# 
a <- ex_range_dpi %>%
  ggplot(aes(x = Period, y = gap, group = Sex, color = Sex)) + 
  geom_line(size = 1.5) + 
  geom_point(size = 4) +
  labs(x = "Period", y = "Range", title = "Range of life expectancy") +
  # labs(x = "Period", y = "Ratio", title = "Ratio of life expectancy") +
  theme_classic() +
  theme(axis.text = element_text(size = 12), axis.title =  element_text(size = 14), legend.text = element_text(size = 12), 
        plot.title = element_text(face = "bold")) +
  theme(axis.text.x = element_text(angle = 50, hjust=1), legend.position = "none")+
  scale_color_manual(labels = c("Female" ,"Male"),
                     values=c("#FF6000", "#339999")) + 
  scale_y_continuous(labels = label_number(accuracy = 0.01)) +
  theme(panel.grid.major = element_line(colour = "lightgrey"))
a

# ------Range and ratio of lifespan variation------# 

# Vector with ages, to be used inside sdv_func 
Age <- c(0:90)

b <- sdv_range_dpi %>%
  ggplot(aes(x = Period, y = gap, group = Sex, color = Sex)) + 
  geom_line(size = 1.5) + 
  geom_point(size = 4) +
  labs(x = "Period", y = "Range", title = "Range of lifespan variation") +
  # labs(x = "Period", y = "Ratio", title = "Ratio of lifespan variation") +
  theme_classic() +
  theme(axis.text = element_text(size = 12), axis.title =  element_text(size = 14), legend.text = element_text(size = 12), 
        plot.title = element_text(face = "bold")) +
  theme(axis.text.x = element_text(angle = 50, hjust=1), legend.position = "none")+
  scale_color_manual(labels = c("Female" ,"Male"),
                     values=c("#FF6000", "#339999")) + 
  scale_y_continuous(labels = label_number(accuracy = 0.01)) +
  theme(panel.grid.major = element_line(colour = "lightgrey"))
b


# ------Slope/Relative index of inequality for life expectancy------# 

c <- sii_dpi %>%
  filter(measure == "sii") %>%
  ggplot(aes(x = Period, y = value, group = Sex, color = Sex)) + 
  geom_line(size = 1.5) + 
  geom_point(size = 4) +
  labs(x = "Period", y = "SII", title = "SII for age-standarised mortality rates") +
  # labs(x = "Period", y = "RII", title = "RII for age-standarised mortality rates") +
  theme_classic() +
  theme(axis.text = element_text(size = 12), axis.title =  element_text(size = 14), legend.text = element_text(size = 12), 
        plot.title = element_text(face = "bold"), legend.title = element_text(size = 13)) +
  theme(axis.text.x = element_text(angle = 50, hjust=1), legend.position = "none")+
  scale_color_manual(labels = c("Female" ,"Male"),
                     values=c("#FF6000", "#339999")) +
  theme(panel.grid.major = element_line(colour = "lightgrey"))
c


# ------Population Total Variation (PTV)------# 
d <- ptv_dpi %>%
  ggplot(aes(x = Period, y = NPTV, group = Sex, color = Sex)) +
  geom_line(size = 1.5) + 
  geom_point(size = 4) +
  labs(x = "Period", y = "PTV", title = "Population total variation at birth") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12), axis.title =  element_text(size = 13), legend.text = element_text(size = 13), 
        plot.title = element_text(face = "bold"), legend.title = element_text(size = 13),
        axis.text.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(angle = 50, hjust=1), legend.position = "none")+
  scale_color_manual(labels = c("Female" ,"Male"),
                     values=c("#FF6000", "#339999")) + 
  theme(panel.grid.major = element_line(colour = "lightgrey")) 
d

# Print all plots together
ggarrange(a, b, c, d, ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom")

# ---------------------------------------------------------------------------- #
#     3. Save results
# ---------------------------------------------------------------------------- #
# ggsave("Graphs/ENG_absolute.pdf", width = 8.98, height = 6.54)

