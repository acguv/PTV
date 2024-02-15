
# ---------------------------------------------------------------------------- #
# Title:   Population total variation
# Country: Motivation figures
# Author:  Gomez-Ugarte Ana C 
# ---------------------------------------------------------------------------- #

# Content:
#   0. Working directory, package and functions
#   1. Read the data
#   2. Preparing the data
#   3. Main figures

# ---------------------------------------------------------------------------- #
#     0. Working directory, package and functions
# ---------------------------------------------------------------------------- #

# To clear everything in R, before start the analysis
rm(list = ls()) 

# Libraries
library(colorspace)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(scales)

# Working directory from safe environment
#setwd("(...)")

# Import functions
source("R/Functions_PTV.R")

# ---------------------------------------------------------------------------- #
#     1. Read the data
# ---------------------------------------------------------------------------- #

# File with fictional exposures
E <- read.delim("Data/E.txt", header = TRUE, sep = "")
E <- as.vector(unlist(E))

# WHO standard population by single age 0-110+
who_std <- read.csv("Data/WHO_std_single_age_110.csv")

# ---------------------------------------------------------------------------- #
#     2. Preparing the data
# ---------------------------------------------------------------------------- #

x <- 0:110
m <- length(x)

# ------Data for motivation figure for 2 groups------# 

# Defining the parameters for each curve
a1 <- 0.0008
b1 <- 0.09
f1 <- GompFX(c(a1,b1),x)

a2 <- 0.00013
b2 <- 0.09
f2 <- GompFX(c(a2,b2),x)

a3 <- 0.0019
b3 <- 0.05
f3 <- GompFX(c(a3,b3),x)

a4 <- 0.000654
b4 <- 0.05
f4 <- GompFX(c(a4,b4),x)

# Estimate life tables, life expectancy and lifespan variation
LT1 <- lifetable.mx(x=seq(0.5,110.5),mx=GompMU(pars=c(a1,b1),ages=x))
LT2 <- lifetable.mx(x=seq(0.5,110.5),mx=GompMU(pars=c(a2,b2),ages=x))
LT3 <- lifetable.mx(x=seq(0.5,110.5),mx=GompMU(pars=c(a3,b3),ages=x))
LT4 <- lifetable.mx(x=seq(0.5,110.5),mx=GompMU(pars=c(a4,b4),ages=x))

e01 <- LT1$ex[1]; e01
e02 <- LT2$ex[2]; e02
e03 <- LT3$ex[1]; e03
e04 <- LT4$ex[2]; e04

sdv1 <- sdv_func_age(LT1)[1]; sdv1
sdv2 <- sdv_func_age(LT2)[1]; sdv2
sdv3 <- sdv_func_age(LT3)[1]; sdv3
sdv4 <- sdv_func_age(LT4)[1]; sdv4

# Estimating the total variation
tvd.12 <- tvd(LT1$dx, LT2$dx)/200000
tvd.34 <- tvd(LT3$dx, LT4$dx)/200000

# Extending the exposures until age 110
e <- c(E, 9,7,6,5,4,3,2,2,1,1,1)

# Estimating the mortality rate of the total population (we assume all groups have equal sizes)
LT1$Dx <- e * LT1$mx
LT2$Dx <- e * LT2$mx
LT3$Dx <- e * LT3$mx
LT4$Dx <- e * LT4$mx

mx_tot.12 <-  (LT1$Dx + LT2$Dx)/(2*e)
LT_tot.12 <- lifetable.mx(x=seq(0.5,110.5),mx=mx_tot.12)
LT_tot.12$dx[LT_tot.12$dx<0] <- 0

mx_tot.34 <-  (LT3$Dx + LT4$Dx)/(2*e)
LT_tot.34 <- lifetable.mx(x=seq(0.5,110.5),mx=mx_tot.34)

# Creating data frame with the age-at-death distributions
f_2groups <- as.data.frame(cbind(f1,f2,f3,f4, tot.12 = LT_tot.12$dx/100000, tot.34 = LT_tot.34$dx/100000))

# PTV 
PTV.12 <- tvd(LT1$dx, LT_tot.12$dx)/200000*.5 + tvd(LT2$dx, LT_tot.12$dx)/200000*.5
PTV.34 <- tvd(LT3$dx, LT_tot.34$dx)/200000*.5 + tvd(LT4$dx, LT_tot.34$dx)/200000*.5


# ------Data for motivation figure for 3 groups------# 

# Defining the parameters for each curve
a5 <- 0.0025
b5 <- 0.09
f5 <- GompFX(c(a5,b5),x)

a6 <- 0.0005
b6 <- 0.09
f6 <- GompFX(c(a6,b6),x)

a7 <- 0.000125
b7 <- 0.09
f7 <- GompFX(c(a7,b7),x)

a8 <- 0.00015
b8 <- 0.09 
f8 <- GompFX(c(a8,b8),x)

# Estimate life tables, life expectancy and lifespan variation
LT5 <- lifetable.mx(x=seq(0.5,110.5),mx=GompMU(pars=c(a5,b5),ages=x))
LT6 <- lifetable.mx(x=seq(0.5,110.5),mx=GompMU(pars=c(a6,b6),ages=x))
LT7 <- lifetable.mx(x=seq(0.5,110.5),mx=GompMU(pars=c(a7,b7),ages=x))
LT8 <- lifetable.mx(x=seq(0.5,110.5),mx=GompMU(pars=c(a8,b8),ages=x))

sdv5 <- sdv_func_age(LT5)[1]; sdv5
sdv6 <- sdv_func_age(LT6)[1]; sdv6
sdv7 <- sdv_func_age(LT7)[1]; sdv7
sdv8 <- sdv_func_age(LT8)[1]; sdv8

# Estimating the mortality rate of the total population (we assume all groups have equal sizes)
LT5$Dx <- e* LT5$mx
LT6$Dx <- e* LT6$mx
LT7$Dx <- e* LT7$mx
LT8$Dx <- e* LT8$mx

mx_tot.57 <-  (LT5$Dx + LT6$Dx + LT7$Dx)/(3*e)
LT_tot.57 <- lifetable.mx(x=seq(0.5,110.5), mx=mx_tot.57)
LT_tot.57$dx[LT_tot.57$dx<0] <- 0

mx_tot.58 <-  (LT5$Dx + LT8$Dx + LT7$Dx)/(3*e)
LT_tot.58 <- lifetable.mx(x=seq(0.5,110.5), mx=mx_tot.58)
LT_tot.58$dx[LT_tot.58$dx<0] <- 0

f_3groups <- as.data.frame(cbind(f5, f6, f7, f8, tot.57 = LT_tot.57$dx/100000, tot.58 = LT_tot.58$dx/100000))

# PTV 
PTV.57 <- tvd(LT5$dx, LT_tot.57$dx)/200000*1/3 + tvd(LT6$dx, LT_tot.57$dx)/200000*1/3 + tvd(LT7$dx, LT_tot.57$dx)/200000*1/3
PTV.58 <- tvd(LT5$dx, LT_tot.58$dx)/200000*1/3 + tvd(LT8$dx, LT_tot.58$dx)/200000*1/3 + tvd(LT7$dx, LT_tot.58$dx)/200000*1/3

# Slope index of inequality of the age-standardized mortality rates
f_who <- cbind(c(0:110), f_3groups[,-c(5,6)])
colnames(f_who)[1] <- "Age"

f_who <- merge(f_who, who_std, by = "Age")

asmr <- f_who %>%  
  gather("f", "mx", f5:f8) %>%
  mutate(mx_w = mx*Prop) %>%
  group_by(f) %>% 
  mutate(asmr = sum(mx_w)*100000)  %>%
  select(f, asmr) %>%
  unique() 

reg_asmr.57 <- lm(c(asmr$asmr[-4]) ~ c(1/6,1/2,5/6))
sii_asmr.57 <- -reg_asmr.57$coefficients[2]
rii_asmr.57 <- (-reg_asmr.57$coefficients[1]-reg_asmr.57$coefficients[2])/(-reg_asmr.57$coefficients[1])

reg_asmr.58 <- lm(c(asmr$asmr[c(1,4,3)]) ~ c(1/6,1/2,5/6))
sii_asmr.58 <- -reg_asmr.58$coefficients[2]
rii_asmr.58 <- (-reg_asmr.58$coefficients[1]-reg_asmr.58$coefficients[2])/(-reg_asmr.58$coefficients[1])


# ------Simulating different weights for the 2 groups example------# 
j = 1
ptv = 0
ptv_w = 0
for (i in seq(0,1,.01)) {
  LT1$Dx <- e * i * LT1$mx
  LT2$Dx <- e * (1-i) * LT2$mx
  mx_tot.12 <-  (LT1$Dx + LT2$Dx)/(e)
  LT_tot.12 <- lifetable.mx(x=seq(0.5,110.5),mx=mx_tot.12)
  temp <- tvd(LT1$dx, LT_tot.12$dx)/200000 + tvd(LT2$dx, LT_tot.12$dx)/200000
  temp_w <- tvd(LT1$dx, LT_tot.12$dx)/200000*i + tvd(LT2$dx, LT_tot.12$dx)/200000*(1-i)
  ptv <- c(ptv, temp)
  ptv_w <- c(ptv_w, temp_w)
  j = j+1
}

wts_sim <- cbind.data.frame(pop_share = rep(seq(0,1,.01),2), method = rep(c("Not weighted", "Weighted"), each = 101), value = c(ptv = ptv[-1], ptv_w = ptv_w[-1]))

# ------Simulating different weights for the 3 groups example------# 
j = 1
ptv3 = 0
ptv3_w = 0
ks = 0
is = 0
for (k in seq(0,1,by=0.05)) {
  for (i in seq(0,1,.01)) {
    if(k+i > .99){
      next
    }
    LT5$Dx <- e * i * LT5$mx
    LT6$Dx <- e * k * LT6$mx
    LT7$Dx <- e * (1-i-k) * LT7$mx    
    mx_tot <-  (LT5$Dx + LT6$Dx + LT7$Dx)/(e)
    LT_tot <- lifetable.mx(x=seq(0.5,110.5), mx=mx_tot)
    temp3 <- tvd(LT5$dx, LT_tot$dx)/200000 + tvd(LT6$dx, LT_tot$dx)/200000 + tvd(LT7$dx, LT_tot$dx)/200000
    temp3_w <-tvd(LT5$dx, LT_tot$dx)/200000*i + tvd(LT6$dx, LT_tot$dx)/200000*k + tvd(LT7$dx, LT_tot$dx)/200000*(1-i-k)
    ptv3 <- c(ptv3, temp3)
    ptv3_w <- c(ptv3_w, temp3_w)
    ks <- c(ks, k)
    is <- c(is, i)
    j = j+1
  }
}

wts3_sim <- cbind.data.frame(is = is[-1], ks = ks[-1], js = (1-is[-1]-ks[-1]),
                             method = rep(c("Not weighted", "Weighted"), each = 1045), 
                             value = c(ptv = ptv3[-1], ptv_w = ptv3_w[-1]))

# ---------------------------------------------------------------------------- #
#     3. Main figures
# ---------------------------------------------------------------------------- #

# ------Panel A Motivation figure------# 
a <- f_2groups %>% 
  ggplot(aes(x = x, y =f1)) +
  geom_ribbon(aes(ymax = f1, ymin = 0), fill="#CC6677", alpha=.5) +
  geom_ribbon(aes(ymax = f2, ymin = 0), fill="#44AA99", alpha=.5) +
  geom_line(aes(x = x, y = f1), size = 1.5, color = "#CC6677") +
  geom_line(aes(x = x, y = f2), size = 1.5, color = "#44AA99") +
  geom_line(aes(x = x, y = tot.12), size = 1, color = "black", linetype = 5) +
  theme_minimal() +
  ggtitle("A") + xlab("") + ylab("Life table deaths (dx)") +
  ylim(c(0, 0.035)) + 
  theme(axis.text = element_text(size = 14), axis.title =  element_text(size = 14)) +  
  annotate("text", x = 4, y = 0.018, label = paste("e[0] ==", sprintf('%.1f', e01)) , size = 4, parse = TRUE, hjust = 0) +
  annotate("text", x = 4, y = 0.015, label = paste("sigma ==", sprintf('%.1f', sdv1)) , size = 4, parse = TRUE, hjust = 0) +
  annotate("text", x = 19, y = 0.031, label = paste("e[0] ==", sprintf('%.1f', e02)) , size = 4, parse = TRUE, hjust = 0) +
  annotate("text", x = 19, y = 0.028, label = paste("sigma ==", sprintf('%.1f', sdv2)) , size = 4, parse = TRUE, hjust = 0) +
  # annotate("text", x = 85, y = 0.035, label = paste("Delta~e[", 0, "] ==",sprintf('%.1f', e02-e01)) , size = 5, parse = TRUE, hjust = 0) +
  annotate("text", x = 85, y = 0.035, label = paste("ratio~e[", 0, "] ==",sprintf('%.1f', e02/e01)) , size = 5, parse = TRUE, hjust = 0) +
  # annotate("text", x = 85, y = 0.031, label = paste("Delta~sigma ==",sprintf('%.1f', sdv2-sdv1)) , size = 5, parse = TRUE, hjust = 0) +
  annotate("text", x = 85, y = 0.031, label = paste("ratio~sigma ==",sprintf('%.1f', sdv2/sdv1)) , size = 5, parse = TRUE, hjust = 0) +
  annotate("text", x = 85, y = 0.026, label = paste("PTV =", sprintf('%.2f',PTV.12)) , size = 5, hjust = 0) +
  geom_segment(aes(x = 42, y = .014, xend = 20, yend = .016),
               arrow = arrow(length = unit(0.3, "cm")), color = "#CC6677", lwd = 1) +
  geom_segment(aes(x = 75, y = .027, xend = 38, yend = .03),
               arrow = arrow(length = unit(0.3, "cm")), color = "#44AA99", lwd = 1) +
  scale_x_continuous(breaks = seq(0, 100, 20))
a

# ------Panel B Motivation figure------# 
b <- f_2groups %>% 
  ggplot(aes(x = x, y = f3)) +
  geom_ribbon(aes(ymax = f3, ymin = 0), fill="#CC6677", alpha=.5) +
  geom_ribbon(aes(ymax = f4, ymin = 0), fill="#44AA99", alpha=.5) +
  geom_line(aes(x = x, y = f3), size = 1.5, color = "#CC6677") +
  geom_line(aes(x = x, y = f4), size = 1.5, color = "#44AA99") +
  geom_line(aes(x = x, y = tot.34), size = 1, color = "black", linetype = 5) +
  theme_minimal() +
  theme(axis.text = element_text(size = 14), axis.title =  element_text(size = 14)) +  
  ggtitle("B") + xlab("") + ylab("") +
  ylim(c(0, 0.035)) + 
  annotate("text", x = 6, y = 0.013, label = paste("e[0] ==", sprintf('%.1f', e03)) , size = 4, parse = TRUE, hjust = 0) +
  annotate("text", x = 6, y = 0.01, label = paste("sigma ==", sprintf('%.1f', sdv3)) , size = 4, parse = TRUE, hjust = 0) +
  annotate("text", x = 40, y = 0.026, label = paste("e[0] ==", sprintf('%.1f', e04)) , size = 4, parse = TRUE, hjust = 0) +
  annotate("text", x = 40, y = 0.023, label = paste("sigma ==", sprintf('%.1f', sdv4)) , size = 4, parse = TRUE, hjust = 0) +
  # annotate("text", x = 85, y = 0.035, label = paste("Delta~e[", 0, "] ==",sprintf('%.1f', e04-e03)) , size = 5, parse = TRUE, hjust = 0) +
  annotate("text", x = 85, y = 0.035, label = paste("ratio~e[", 0, "] ==",sprintf('%.1f', e04/e03)) , size = 5, parse = TRUE, hjust = 0) +
  # annotate("text", x = 85, y = 0.031, label = paste("Delta~sigma ==",sprintf('%.1f', sdv4-sdv3)) , size = 5, parse = TRUE, hjust = 0) +
  annotate("text", x = 85, y = 0.031, label = paste("ratio~sigma ==",sprintf('%.1f', sdv4/sdv3)) , size = 5, parse = TRUE, hjust = 0) +
  annotate("text", x = 85, y = 0.026, label = paste("PTV =", sprintf('%.2f', PTV.34)) , size = 5, hjust = 0) +
  geom_segment(aes(x = 45, y = .01, xend = 24, yend = .011),
               arrow = arrow(length = unit(0.3, "cm")), color = "#CC6677", lwd = 1) +
  geom_segment(aes(x = 88, y = .015, xend = 56, yend = .025),
               arrow = arrow(length = unit(0.3, "cm")), color = "#44AA99", lwd = 1) +
  scale_x_continuous(breaks = seq(0, 100, 20))
b

# ------Panel c Motivation figure------# 
c <- f_3groups %>% 
  ggplot(aes(x = x, y = f5)) +
  geom_ribbon(aes(ymax = f5, ymin = 0), fill="#CC6677", alpha=.5) +
  geom_ribbon(aes(ymax = f6, ymin = 0), fill="#44AA99", alpha=.5) +
  geom_ribbon(aes(ymax = f7, ymin = 0), fill="#EE7733", alpha=.5) +
  geom_line(aes(x = x, y = f5), size = 1.5, color = "#CC6677") + 
  geom_line(aes(x = x, y = f6), size = 1.5, color = "#44AA99") + 
  geom_line(aes(x = x, y = f7), size = 1.5, color = "#EE7733") + 
  geom_line(aes(x = x, y = tot.57), size = 1, color = "black", linetype = 5) + 
  theme_minimal() +
  theme(axis.text = element_text(size = 14), axis.title =  element_text(size = 14)) +
  ylim(c(0,0.035)) + 
  ggtitle("C") + xlab("Age") + ylab("Life table deaths (dx)") +
  # annotate("text", x = 85, y = 0.034, label = paste("SII =", sprintf('%.1f', sii_asmr.57)) , size = 5, hjust = 0) +
  annotate("text", x = 85, y = 0.034, label = paste("RII =", sprintf('%.2f', rii_asmr.57)) , size = 5, hjust = 0) +
  annotate("text", x = 85, y = 0.030, label = paste("PTV =", sprintf('%.2f', PTV.57)) , size = 5, hjust = 0) +
  scale_x_continuous(breaks = seq(0, 100, 20))
c

# ------Panel D Motivation figure------# 
d <- f_3groups %>% 
  ggplot(aes(x = x, y = f5)) +
  geom_ribbon(aes(ymax = f5, ymin = 0), fill="#CC6677", alpha=.5) +
  geom_ribbon(aes(ymax = f8, ymin = 0), fill="#44AA99", alpha=.5) +
  geom_ribbon(aes(ymax = f7, ymin = 0), fill="#EE7733", alpha=.5) +
  geom_line(aes(x = x, y = f5), size = 1.5, color = "#CC6677") + 
  geom_line(aes(x = x, y = f8), size = 1.5, color = "#44AA99") + 
  geom_line(aes(x = x, y = f7), size = 1.5, color = "#EE7733") + 
  geom_line(aes(x = x, y = tot.58), size = 1, color = "black", linetype = 5) + 
  theme_minimal() +
  theme(axis.text = element_text(size = 14), axis.title =  element_text(size = 14)) +
  ylim(c(0,0.035)) + 
  ggtitle("D") + xlab("Age") + ylab("") +
  # annotate("text", x = 85, y = 0.034, label = paste("SII =", sprintf('%.1f', sii_asmr.58)) , size = 5, hjust = 0) +
  annotate("text", x = 85, y = 0.034, label = paste("RII =", sprintf('%.2f', rii_asmr.58)) , size = 5, hjust = 0) +
  annotate("text", x = 85, y = 0.030, label = paste("PTV =", sprintf('%.2f', PTV.58)) , size = 5, hjust = 0) +
  scale_x_continuous(breaks = seq(0, 100, 20))
d

# ------Motivation figure------# 
ggarrange(a,b,c,d, ncol = 2, nrow = 2)
# ggsave("Graphs/Motivation.pdf", width = 10.35, height = 6.91)


# ------Supplementary Material figures of the effect of the weights on the PTV for 2 groups------# 
ggplot(wts_sim[wts_sim$method == "Weighted",], aes(x = pop_share, y = value)) +
  geom_line(lwd = 1.5) +
  # facet_grid(~ method, scales="free") + 
  geom_point(data = data.frame(pop_share = 0, value = ptv_w[1], method = "Weighted"), size = 4) +
  geom_point(data = data.frame(pop_share = 1, value = ptv_w[102], method = "Weighted"), size = 4) +
  xlab("Share of population in the worst-off group") + ylab("PTV") + 
  theme_bw() +
  scale_x_continuous(labels = label_number(accuracy = 0.1)) +
  theme(strip.text = element_text(size = 12), axis.title = element_text(size = 11), axis.text =  element_text(size = 10))

# ggsave("Graphs/PTV_weights_2groups.pdf", width = 5.02, height = 3.23)


# ------Supplementary Material figures of the effect of the weights on the PTV for 3 groups------# 
ggplot(wts3_sim[wts3_sim$method == "Weighted",], aes(x = is, y = value, color = as.factor(ks))) +
  geom_line(lwd = 1) +
  # facet_grid(~ method, scales="free") + 
  xlab("Share of population in the worst-off group") + ylab("PTV") + 
  theme_bw() +
  scale_x_continuous(labels = label_number(accuracy = 0.1)) +
  theme(strip.text = element_text(size = 12), axis.title = element_text(size = 11), axis.text =  element_text(size = 10)) + 
  scale_color_discrete(name = "Weights of\nmiddle group")+ 
  guides(color=guide_legend(ncol=2))

# ggsave("Graphs/PTV_weights_3groups.pdf", width = 7.74, height = 4.36)





