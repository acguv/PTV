
# ---------------------------------------------------------------------------- #
# Paper:   Reassessing socioeconomic inequalities in mortality via distributional similarities
# Title:   Plot results
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
load("Results/ENG_indicators.RData")

# 2. Standard population
who_std <- read.csv("Data/WHO_std_single_age_90+.csv")
names(who_std)[1] <- "Age"

# ---------------------------------------------------------------------------- #
#     2. Main figures
# ---------------------------------------------------------------------------- #

# ------Range and ratio of life expectancy------# 
a <- ex_range_dpi %>%
  ggplot(aes(x = Period, y = ratio, group = Sex, color = Sex)) + 
  geom_line(size = 1.5) + 
  geom_point(size = 4) +
  # labs(x = "Period", y = "Range", title = "Range of life expectancy") +
  labs(x = "Period", y = "Ratio", title = "Ratio of life expectancy") +
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
  ggplot(aes(x = Period, y = ratio, group = Sex, color = Sex)) + 
  geom_line(size = 1.5) + 
  geom_point(size = 4) +
  # labs(x = "Period", y = "Range", title = "Range of lifespan variation") +
  labs(x = "Period", y = "Ratio", title = "Ratio of lifespan variation") +
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
  filter(measure == "rii") %>%
  ggplot(aes(x = Period, y = value, group = Sex, color = Sex)) + 
  geom_line(size = 1.5) + 
  geom_point(size = 4) +
  # labs(x = "Period", y = "SII", title = "SII for age-standarised mortality rates") +
  labs(x = "Period", y = "RII", title = "RII for age-standarised mortality rates") +
  theme_classic() +
  theme(axis.text = element_text(size = 12), axis.title =  element_text(size = 14), legend.text = element_text(size = 12), 
        plot.title = element_text(face = "bold"), legend.title = element_text(size = 13)) +
  theme(axis.text.x = element_text(angle = 50, hjust=1), legend.position = "none")+
  scale_color_manual(labels = c("Female" ,"Male"),
                     values=c("#FF6000", "#339999")) +
  theme(panel.grid.major = element_line(colour = "lightgrey"))
c

# ------Pairwise non-overlap index (NOI)------# 
d <- noi_pair_dpi %>%
  ggplot(aes(x = Period, y = noi_pair, group = Sex, color = Sex)) +
  geom_line(size = 1.5) + 
  geom_point(size = 4) +
  labs(x = "Period", title = "Pairwise non-overlap index at birth") +
  ylab(expression(S^P)) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12), axis.title =  element_text(size = 13), legend.text = element_text(size = 13), 
        plot.title = element_text(face = "bold"), legend.title = element_text(size = 13),
        axis.text.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(angle = 50, hjust=1), legend.position = "none")+
  scale_color_manual(labels = c("Female" ,"Male"),
                     values=c("#FF6000", "#339999")) + 
  theme(panel.grid.major = element_line(colour = "lightgrey")) 
d

# ------Total non-overlap index (NOI)------# 
e <- noi_t_dpi %>%
  ggplot(aes(x = Period, y = noi_t, group = Sex, color = Sex)) +
  geom_line(size = 1.5) + 
  geom_point(size = 4) +
  labs(x = "Period", title = "Total non-overlap index at birth") +
  ylab(expression(S^T)) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12), axis.title =  element_text(size = 13), legend.text = element_text(size = 13), 
        plot.title = element_text(face = "bold"), legend.title = element_text(size = 13),
        axis.text.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(angle = 50, hjust=1), legend.position = "none")+
  scale_color_manual(labels = c("Female" ,"Male"),
                     values=c("#FF6000", "#339999")) + 
  theme(panel.grid.major = element_line(colour = "lightgrey")) 
e

# ------Stratification index (Zhou et al.)------# 
f <- S_dpi %>%
  ggplot(aes(x = Period, y = S, group = Sex, color = Sex)) +
  geom_line(size = 1.5) + 
  geom_point(size = 4) +
  labs(x = "Period", title = "Stratification index at birth", y = "S") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12), axis.title =  element_text(size = 13), legend.text = element_text(size = 13), 
        plot.title = element_text(face = "bold"), legend.title = element_text(size = 13),
        axis.text.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(angle = 50, hjust=1), legend.position = "none")+
  scale_color_manual(labels = c("Female" ,"Male"),
                     values=c("#FF6000", "#339999")) + 
  theme(panel.grid.major = element_line(colour = "lightgrey")) 
f


# ------Pairwise outsurvival probability------# 
g <- OV_pair_dpi %>%
  ggplot(aes(x = Period, y = OV_pair, group = Sex, color = Sex)) +
  geom_line(size = 1.5) + 
  geom_point(size = 4) +
  labs(x = "Period", title = "Pairwise outsurvival probability at birth", y = "OV") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12), axis.title =  element_text(size = 13), legend.text = element_text(size = 13), 
        plot.title = element_text(face = "bold"), legend.title = element_text(size = 13),
        axis.text.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(angle = 50, hjust=1), legend.position = "none")+
  scale_color_manual(labels = c("Female" ,"Male"),
                     values=c("#FF6000", "#339999")) + 
  theme(panel.grid.major = element_line(colour = "lightgrey")) 
g


# ------Population Total Variation (PTV) (for working paper)------# 
h <- ptv_dpi %>%
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
h

# Print all plots together
ggarrange(a, b, c, d, ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom")

# ---------------------------------------------------------------------------- #
#     3. Save results
# ---------------------------------------------------------------------------- #
# ggsave("Graphs/Fig5.pdf", width = 8.98, height = 6.54)

