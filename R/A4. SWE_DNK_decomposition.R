
# ---------------------------------------------------------------------------- #
# Paper:   Reassessing socioeconomic inequalities in mortality via distributional similarities
# Title:   Stepwise decomposition of inequality measures
# Country: Denmark and Sweden
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
#     1. Read data
# ---------------------------------------------------------------------------- #

# Read data
load("Data/SWE_DNK_LT_tot.RData")
load("Data/SWE_DNK_Edu_ranks.RData")
load("Data/SWE_DNK_Edu_weigths.RData")

who_std <- read.csv("Data/WHO_std_5_age.csv")

# ---------------------------------------------------------------------------- #
#     2. Decomposition
# ---------------------------------------------------------------------------- #

# Function to get the the mortality rates in each group for a specific sex, country and period
get_mx = function(sex, period, country){
  mx <- lt_edu %>%  
    filter(Education != "Total") %>%
    # We don't need the mortality rates of the total population because we construct it from the other
    select(Country, Period, Sex, AgeGrp, Age, Education, mx) %>%
    arrange(Country, Sex, Period, Education) %>% 
    filter(Sex == sex, Period == period & Country == country) %>%
    .$mx
  return(mx)
}

# Function to get the populations shares in each group for a specific sex, country and period
get_weights = function(sex, period, country){
  wts <- lt_edu %>%
    filter(Education != "Total") %>%
    arrange(Country, Sex, Period, Education) %>%
    filter(Sex == sex, Period == period & Country == country) %>%
    .$EduPop
  return(wts)
}

# Function to order the decomposition results
contr_fun <- function(decomp_results, country, sex){
  contr <- data.frame(cbind(matrix(decomp_results[1:39], ncol = 3), 
                            rowSums(matrix(decomp_results[40:78], ncol = 3))))
  names(contr) <- c("Low", "Middle", "High", "Pop. Composition")
  rel_contr <- contr/sum(contr)
  names(rel_contr) <- c("Low", "Middle", "High", "Pop. Composition")
  contr <- rbind(contr, rel_contr)
  contr %>%
    mutate(Country = country, 
           Sex = sex,
           Series = rep(c("Absolute", "Relative"), each = 13),
           Age = rep(seq(30,90,5),2)) %>%
    gather("Education", "Contribution", "Low":"Pop. Composition") 
}

# Get weight for each country and sex
wts_den_m_91 <- get_weights(sex = "Male", period = "1991-1995", country = "Denmark")
wts_den_m_11 <- get_weights(sex = "Male", period = "2011-2015", country = "Denmark")

wts_den_f_91 <- get_weights(sex = "Female", period = "1991-1995", country = "Denmark")
wts_den_f_11 <- get_weights(sex = "Female", period = "2011-2015", country = "Denmark")

wts_swe_m_91 <- get_weights(sex = "Male", period = "1991-1995", country = "Sweden")
wts_swe_m_11 <- get_weights(sex = "Male", period = "2011-2015", country = "Sweden")

wts_swe_f_91 <- get_weights(sex = "Female", period = "1991-1995", country = "Sweden")
wts_swe_f_11 <- get_weights(sex = "Female", period = "2011-2015", country = "Sweden")


# Get dx for each country and sex
mx_den_m_91 <- get_mx(sex = "Male", period = "1991-1995", country = "Denmark")
mx_den_m_11 <- get_mx(sex = "Male", period = "2011-2015", country = "Denmark")

mx_den_f_91 <- get_mx(sex = "Female", period = "1991-1995", country = "Denmark")
mx_den_f_11 <- get_mx(sex = "Female", period = "2011-2015", country = "Denmark")

mx_swe_m_91 <- get_mx(sex = "Male", period = "1991-1995", country = "Sweden")
mx_swe_m_11 <- get_mx(sex = "Male", period = "2011-2015", country = "Sweden")

mx_swe_f_91 <- get_mx(sex = "Female", period = "1991-1995", country = "Sweden")
mx_swe_f_11 <- get_mx(sex = "Female", period = "2011-2015", country = "Sweden")


# Create vectors for the decomposition
mx_w_den_m_91 <- c(mx_den_m_91, wts_den_m_91)
mx_w_den_m_11 <- c(mx_den_m_11, wts_den_m_11)

mx_w_den_f_91 <- c(mx_den_f_91, wts_den_f_91)
mx_w_den_f_11 <- c(mx_den_f_11, wts_den_f_11)

mx_w_swe_m_91 <- c(mx_swe_m_91, wts_swe_m_91)
mx_w_swe_m_11 <- c(mx_swe_m_11, wts_swe_m_11)

mx_w_swe_f_91 <- c(mx_swe_f_91, wts_swe_f_91)
mx_w_swe_f_11 <- c(mx_swe_f_11, wts_swe_f_11)

#-----------Decomposition of the PAIRWISE NOI---------# 
# Function to estimate the pairwise NOI given vectors with the education-specific mortality rates (mx)
# and the population shares in each subgroup (w)
noi_pair_func <- function(mx_w, n = 3){
  w <- matrix(as.matrix(mx_w[(length(mx_w)/2+1):length(mx_w)]), ncol = n)
  w_edu <- colSums(w)/sum(w)
  mx <- matrix(as.matrix(mx_w[1:(length(mx_w)/2)]), ncol = n)
  dx <- apply(mx, 2, function(x) lifetable.mx(seq(30,90,5), x, ax=NULL)$dx)
  # dx <- cbind(dx, dx_tot)
  dx.d <- dx/colSums(dx)
  comb_dx <- combn(ncol(dx),2)
  
  Sij <- rep(NA,ncol(comb_dx))
  weights <- rep(NA,ncol(comb_dx))
  i <- 1
  for (i in 1:ncol(comb_dx)){
    Sij[i] <- noi_func(dx[,comb_dx[1,i]],dx[,comb_dx[2,i]])
    weights[i] <- w_edu[comb_dx[1,i]]*w_edu[comb_dx[2,i]]
  }
  
  noi_t <- sum(Sij*weights)/sum(weights)
  noi_t
}

# Decomposition
decomp_den_m <- stepwise_replacement(noi_pair_func, mx_w_den_m_91, mx_w_den_m_11) 
decomp_den_f <- stepwise_replacement(noi_pair_func, mx_w_den_f_91, mx_w_den_f_11) 

decomp_swe_m <- stepwise_replacement(noi_pair_func, mx_w_swe_m_91, mx_w_swe_m_11) 
decomp_swe_f <- stepwise_replacement(noi_pair_func, mx_w_swe_f_91, mx_w_swe_f_11) 


# Ordering the decomposition results
contr_den_m <- contr_fun(decomp_den_m, "Denmark", "Males")
contr_den_f <- contr_fun(decomp_den_f, "Denmark", "Females")
contr_swe_m <- contr_fun(decomp_swe_m, "Sweden", "Males")
contr_swe_f <- contr_fun(decomp_swe_f, "Sweden", "Females")

contr_den_m %>%
  group_by(Series, Education) %>%
  summarise(sum(Contribution))

# Creating single data frame with all the decomposition results
contr_all_noi_pair <- rbind(contr_den_f, contr_den_m, contr_swe_f, contr_swe_m)

# Graph of the decomposition results by age
decomp_noi_pair_graph <- contr_all_noi_pair %>%
  filter(Series == "Absolute") %>%
  mutate(Education = fct_relevel(Education, c("Low", "Middle","High", "Pop. Composition"))) %>%
  ggplot(aes(x= Age, y = Contribution, group = Education, fill = Education)) +
  geom_bar(stat = "identity", color = "black") +
  facet_grid(Country ~ Sex) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12), axis.title =  element_text(size = 13, face = "bold"), legend.text = element_text(size = 13), 
        plot.title = element_text(face = "bold"), legend.title = element_text(size = 13),
        axis.text.y = element_text(size = 12)) + xlab("Age") + ylab(expression(paste("Change in ", S^P))) + 
  scale_fill_viridis_d() +
  theme(panel.grid.major.y = element_line(colour = "lightgrey")) +
  scale_x_continuous(breaks = seq(0, 90, by = 10))

decomp_noi_pair_graph

#------Decomposition of the SII---------#

# Function to estimate the SII given vectors with the education-specific mortality rates (mx) 
# and the population shares in each subgroup (w)
who_std <- who_std[c(8:20),] %>%
  mutate(Prop = Prop/sum(Prop))

sii_func <- function(mx_w, n = 3){
  w <- matrix(as.matrix(mx_w[((length(mx_w)/2)+1):length(mx_w)]), ncol = n)
  w_edu <- colSums(w)/sum(w)
  w_ranks <- c(.5*w_edu[1]+w_edu[2]+w_edu[3],.5*w_edu[2]+w_edu[3],.5*w_edu[3])
  mx <- matrix(as.matrix(mx_w[1:(length(mx_w)/2)]), ncol = n)
  asmr <- colSums(mx*who_std$Prop)*100000
  sii <- lm(asmr ~ w_ranks)$coefficients[2]
  sii
}

# Decomposition
decomp_den_m_sii <- stepwise_replacement(sii_func, mx_w_den_m_91, mx_w_den_m_11) 
decomp_den_f_sii <- stepwise_replacement(sii_func, mx_w_den_f_91, mx_w_den_f_11) 

decomp_swe_m_sii <- stepwise_replacement(sii_func, mx_w_swe_m_91, mx_w_swe_m_11) 
decomp_swe_f_sii <- stepwise_replacement(sii_func, mx_w_swe_f_91, mx_w_swe_f_11) 

# Ordering the results of the decomposition
contr_den_m_sii <- contr_fun(decomp_den_m_sii, "Denmark", "Males")
contr_den_f_sii <- contr_fun(decomp_den_f_sii, "Denmark", "Females")
contr_swe_m_sii <- contr_fun(decomp_swe_m_sii, "Sweden", "Males")
contr_swe_f_sii <- contr_fun(decomp_swe_f_sii, "Sweden", "Females")

contr_all_sii <- rbind(contr_den_f_sii, contr_den_m_sii, contr_swe_f_sii, contr_swe_m_sii) 

# Graph of the decomposition by age
decomp_sii_graph <- contr_all_sii %>%
  filter(Series == "Absolute") %>%
  mutate(Education = fct_relevel(Education, c("Low", "Middle","High", "Pop. Composition"))) %>%
  ggplot(aes(x= Age, y = Contribution, group = Education, fill = Education)) +
  geom_bar(stat = "identity", color = "black") +
  facet_grid(Country ~ Sex) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12), axis.title =  element_text(size = 13, face = "bold"), legend.text = element_text(size = 13), 
        plot.title = element_text(face = "bold"), legend.title = element_text(size = 13),
        axis.text.y = element_text(size = 12)) + xlab("Sex") + ylab("Change in SII") + 
  scale_fill_viridis_d() +
  theme(panel.grid.major.y = element_line(colour = "lightgrey")) +
  scale_x_continuous(breaks = seq(0, 90, by = 10))

#------Decomposition of the RII---------#
# Function to estimate the RII given vectors with the education-specific mortality rates (mx) 
# and the population shares in each subgroup (w)
who_std <- who_std[c(8:20),] %>%
  mutate(Prop = Prop/sum(Prop))

rii_func <- function(mx_w, n = 3){
  w <- matrix(as.matrix(mx_w[((length(mx_w)/2)+1):length(mx_w)]), ncol = n)
  w_edu <- colSums(w)/sum(w)
  w_ranks <- c(.5*w_edu[1]+w_edu[2]+w_edu[3],.5*w_edu[2]+w_edu[3],.5*w_edu[3])
  mx <- matrix(as.matrix(mx_w[1:(length(mx_w)/2)]), ncol = n)
  asmr <- colSums(mx*who_std$Prop)*100000
  fit <- lm(asmr ~ w_ranks)
  rii <- (fit$coefficients[1]+fit$coefficients[2])/fit$coefficients[1]
  rii
}

# Decomposition
decomp_den_m_rii <- stepwise_replacement(rii_func, mx_w_den_m_91, mx_w_den_m_11) 
decomp_den_f_rii <- stepwise_replacement(rii_func, mx_w_den_f_91, mx_w_den_f_11) 

decomp_swe_m_rii <- stepwise_replacement(rii_func, mx_w_swe_m_91, mx_w_swe_m_11) 
decomp_swe_f_rii <- stepwise_replacement(rii_func, mx_w_swe_f_91, mx_w_swe_f_11) 

# Ordering the results of the decomposition
contr_den_m_rii <- contr_fun(decomp_den_m_rii, "Denmark", "Males")
contr_den_f_rii <- contr_fun(decomp_den_f_rii, "Denmark", "Females")
contr_swe_m_rii <- contr_fun(decomp_swe_m_rii, "Sweden", "Males")
contr_swe_f_rii <- contr_fun(decomp_swe_f_rii, "Sweden", "Females")

contr_all_rii <- rbind(contr_den_f_rii, contr_den_m_rii, contr_swe_f_rii, contr_swe_m_rii) 

# Graph of the decomposition by age
decomp_rii_graph <- contr_all_rii%>%
  filter(Series == "Absolute") %>%
  mutate(Education = fct_relevel(Education, c("Low", "Middle","High", "Pop. Composition"))) %>%
  ggplot(aes(x= Age, y = Contribution, group = Education, fill = Education)) +
  geom_bar(stat = "identity", color = "black") +
  facet_grid(Country ~ Sex) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12), axis.title =  element_text(size = 13, face = "bold"), legend.text = element_text(size = 13), 
        plot.title = element_text(face = "bold"), legend.title = element_text(size = 13),
        axis.text.y = element_text(size = 12)) + xlab("Sex") + ylab("Change in RII") + 
  scale_fill_viridis_d() +
  theme(panel.grid.major.y = element_line(colour = "lightgrey")) +
  scale_x_continuous(breaks = seq(0, 90, by = 10))


#------Decomposition of the PTV (for working paper)---------#
# Function to estimate the PTV given vectors with the education-specific mortality rates (mx)
# and the population shares in each subgroup (w)
ptv_func <- function(mx_w, n = 3){
  w <- matrix(as.matrix(mx_w[(length(mx_w)/2+1):length(mx_w)]), ncol = n)
  w_edu <- colSums(w)/sum(w)
  mx <- matrix(as.matrix(mx_w[1:(length(mx_w)/2)]), ncol = n)
  mx_tot <- rowSums(mx*w)/rowSums(w)
  dx_tot <- lifetable.mx(seq(30,90,5), mx_tot, ax=NULL)$dx
  dx <- apply(mx, 2, function(x) lifetable.mx(seq(30,90,5), x, ax=NULL)$dx)
  dx <- cbind(dx, dx_tot)
  dx.d <- dx/colSums(dx)
  diffs <- apply(dx.d,2,function(x) x-dx.d[,n+1])
  tv <- apply(diffs,2, function(x) sum(abs(x), na.rm = TRUE)/2)
  ptv <- sum(tv[-(n+1)]*w_edu)
  ptv
}

# Decomposition
decomp_den_m <- stepwise_replacement(ptv_func, mx_w_den_m_91, mx_w_den_m_11) 
decomp_den_f <- stepwise_replacement(ptv_func, mx_w_den_f_91, mx_w_den_f_11) 

decomp_swe_m <- stepwise_replacement(ptv_func, mx_w_swe_m_91, mx_w_swe_m_11) 
decomp_swe_f <- stepwise_replacement(ptv_func, mx_w_swe_f_91, mx_w_swe_f_11) 

# Ordering the decomposition results
contr_den_m <- contr_fun(decomp_den_m, "Denmark", "Males")
contr_den_f <- contr_fun(decomp_den_f, "Denmark", "Females")
contr_swe_m <- contr_fun(decomp_swe_m, "Sweden", "Males")
contr_swe_f <- contr_fun(decomp_swe_f, "Sweden", "Females")

# Creating single data frame with all the decomposition results
contr_all_ptv <- rbind(contr_den_f, contr_den_m, contr_swe_f, contr_swe_m)

# Graph of the decomposition results by age
decomp_ptv_graph <- contr_all_ptv %>%
  filter(Series == "Absolute") %>%
  mutate(Education = fct_relevel(Education, c("Low", "Middle","High", "Pop. Composition"))) %>%
  ggplot(aes(x= Age, y = Contribution, group = Education, fill = Education)) +
  geom_bar(stat = "identity", color = "black") +
  facet_grid(Country ~ Sex) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12), axis.title =  element_text(size = 13, face = "bold"), legend.text = element_text(size = 13), 
        plot.title = element_text(face = "bold"), legend.title = element_text(size = 13),
        axis.text.y = element_text(size = 12)) + xlab("Sex") + ylab("Change in PTV") + 
  scale_fill_viridis_d() +
  theme(panel.grid.major.y = element_line(colour = "lightgrey")) +
  scale_x_continuous(breaks = seq(0, 90, by = 10))


# ---------------------------------------------------------------------------- #
#     3. Save results
# ---------------------------------------------------------------------------- #
# save(contr_all_noi_pair, contr_all_rii, file = "Results/SWE_DNK_decomp_results_stepwise.RData")
