
# ---------------------------------------------------------------------------- #
# Title:   Population total variation
# Country: Denmark and Sweden by education level
# Author:  Gomez-Ugarte Ana C 
# ---------------------------------------------------------------------------- #

# Content:
#   0. Working directory, package and functions
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

# Download data from the supplementary information from N?meth et al. (2021) 
# https://doi.org/10.1186/s12963-021-00264-1

# Read data

# Life tables by country, sex and education
load("Data/SWE_DNK_LT_tot.RData")

# Inequality in mortality indicators
load("Results/SWE_DNK_indicators.RData")


# ---------------------------------------------------------------------------- #
#     2. Main figures
# ---------------------------------------------------------------------------- #


# ------Swedish females age at death distribution------# 
lt_edu %>%
  filter(Country  == "Sweden" & Sex == "Female") %>%
  filter(Period == "1991-1995" | Period == "2011-2015") %>%
  select(Age, AgeGrp, dx.d, Period, Education) %>%
  mutate(Education = fct_relevel(Education, 
                                 "Low", "Middle", "High","Total")) %>%
  ggplot(aes(x = Age, y = dx.d, group = interaction(Period, Education), color = Period, fill = interaction(Period, Education))) +
  geom_col(position = "identity", alpha = .6) + facet_grid(~Education) +
  facet_grid(~Education) +
  scale_color_manual(labels = c("1991-1995", "2011-2015"),
                     values=c("transparent", "black")) +
  scale_fill_manual(labels = c("Low", "Middle"  , "High"),
                    values=c("darkgrey","#CC6677","darkgrey" ,"#44AA99","darkgrey","#EE7733","darkgrey","darkblue"), guide = "none") +  
  theme_classic() +
  theme(axis.text = element_text(size = 10), axis.title =  element_text(size = 12), legend.text = element_text(size = 10), 
        plot.title = element_text(face = "bold")) +
  ggtitle("") + xlab("Age") + ylab("Life table deaths (dx)") +
  scale_y_continuous(label=comma, limits = c(0,.45)) +
  theme(panel.grid.major = element_line(colour = "lightgrey")) +
  scale_x_continuous(breaks = seq(30, 90, 20), labels = c(seq(30, 80, 20), "90+"))

# ggsave("Graphs/SWE_females_ADD.pdf", width = 7.21, height = 3.47)


# ------Range and ratio of life expectancy------# 
a <- ex_range_edu %>%
  ggplot(aes(x = Period, y = ratio,  color = Country, lty = Sex, group = interaction(Country, Sex))) + 
  geom_line(lwd = 1.5) + 
  # labs(x = "", y = "Range", title = "Range of life expectancy") +
  labs(x = "", y = "Ratio", title = "Ratio of life expectancy") +
  theme_classic() +
  theme(axis.text = element_text(size = 12), axis.title =  element_text(size = 14), legend.text = element_text(size = 12), 
        plot.title = element_text(face = "bold")) +
  theme(axis.text.x = element_text(angle = 50, hjust=1), legend.position = "none")+
  scale_color_manual(labels = c("Denmark" ,"Sweden"),
                     values=c("#D41159", "#1A85FF"))+ 
  scale_y_continuous(labels = label_number(accuracy = 0.01)) +
  theme(panel.grid.major = element_line(colour = "lightgrey"), panel.grid.minor = element_line(colour = "lightgrey"))
a


# ------Range and ratio of lifespan variation------# 

# Vector with ages, to be used inside sdv_func 
Age <- seq(30,90,5)

b <- sdv_range_edu %>%
  ggplot(aes(x = Period, y = ratio,  color = Country, lty = Sex, group = interaction(Country, Sex))) + 
  geom_line(size = 1.5) + 
  # labs(x = "", y = "Range", title = "Range of lifespan variation") +
  labs(x = "", y = "Ratio", title = "Ratio of lifespan variation") +
  theme_classic() +
  theme(axis.text = element_text(size = 12), axis.title =  element_text(size = 14), legend.text = element_text(size = 12), 
        plot.title = element_text(face = "bold")) +
  theme(axis.text.x = element_text(angle = 50, hjust=1), legend.position = "none")+
  scale_color_manual(labels = c("Denmark" ,"Sweden"),
                     values=c("#D41159", "#1A85FF")) + 
  scale_y_continuous(labels = label_number(accuracy = 0.01)) +
  theme(panel.grid.major = element_line(colour = "lightgrey"), panel.grid.minor = element_line(colour = "lightgrey"))
b


# ------Slope/Relative index of inequality for life expectancy------# 

# With real population shares by education level
c <- sii_edu %>%
  filter(measure == "rii") %>%
  ggplot(aes(x = Period, y = value, color = Country, lty = Sex, group = interaction(Country, Sex))) + 
  geom_line(size = 1.5) + 
  # labs(x = "Year", y = "SII", title = "SII of age-standarised mortality rates") +
  labs(x = "Year", y = "RII", title = "RII of life expectancy") +
  theme_classic() +
  theme(axis.text = element_text(size = 12), axis.title =  element_text(size = 13), legend.text = element_text(size = 13), 
        plot.title = element_text(face = "bold"), legend.title = element_text(size = 13)) +
  theme(axis.text.x = element_text(angle = 50, hjust=1), legend.position = "none")+
  scale_color_manual(labels = c("Denmark" ,"Sweden"),
                     values=c("#D41159", "#1A85FF")) +
  theme(panel.grid.major = element_line(colour = "lightgrey"), panel.grid.minor = element_line(colour = "lightgrey"))
c


# ------Population Total Variation (PTV)------# 
d <- ptv_edu %>%
  ggplot(aes(x = Period, y = NPTV, color = Country, lty = Sex, group = interaction(Country, Sex))) + 
  geom_line(size = 1.5) + 
  labs(x = "Year", y = "PTV", title = "Population total variation at age 30") +
  theme_classic() +
  theme(axis.text = element_text(size = 12), axis.title =  element_text(size = 13), legend.text = element_text(size = 13), 
        plot.title = element_text(face = "bold"), legend.title = element_text(size = 13)) +
  theme(axis.text.x = element_text(angle = 50, hjust=1), legend.position = "none")+
  scale_color_manual(labels = c("Denmark" ,"Sweden"),
                     values=c("#D41159", "#1A85FF")) +
  # scale_color_manual(values = c("#1A85FF")) + 
  scale_y_continuous(labels = label_number(accuracy = 0.01)) +
  theme(panel.grid.major = element_line(colour = "lightgrey"), panel.grid.minor = element_line(colour = "lightgrey"))
d

# Print all plots together
ggarrange(a, b, c, d, ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom")

# ---------------------------------------------------------------------------- #
#     3. Save results
# ---------------------------------------------------------------------------- #

# ggsave("Graphs/SWE_DNK_relative.pdf", width = 8.98, height = 6.54)
