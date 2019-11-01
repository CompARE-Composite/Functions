################################################################
# CASE STUDY T2E
# CompARE paper - 2019
# Jordi Cortes 
################################################################
rm(list=ls())

############################################
# Load functions and libraries
############################################
library(ggplot2)
library(snowfall)
library(copula)
library(numDeriv)

source('Rfunctions_T2E.R')


############################################
# Values from OASIS-6 TRIAL
############################################
p10 <- 0.125   # Probability of observing event in endpoint 1 in control arm
p20 <- 0.037   # Probability of observing event in endpoint 2 in control arm
HR1 <- 0.83    # HR endpoint 1
HR2 <- 0.65    # Minimum HR endpoint 2
beta1 <- 1     # Shape parameter for Weibull distribution endpoint 1
beta2 <- 1     # Shape parameter for Weibull distribution endpoint 2
rho <- c(0.01, 0.15, 0.3, 0.5, 0.7, 0.9, 0.95) # Correlation array


############################################
# Calculate ARE and create dataset
############################################
dataset <- Different_scenarios(rho0=rho,beta1=beta1,beta2=beta2,HR1=HR1,HR2=HR2,
                               p1=p10,p2=p20,case=3,copula='Gumbel',rhoType='Spearman')


############################################
# Plot --> ARE vs. Correlation
############################################
ggplot(dataset,aes(x=rho,y=ARE,col=as.factor(HR2))) + geom_line(size=1.3) + geom_hline(yintercept = 1,col='black',linetype=2) +
  scale_y_log10(breaks=c(0.1,0.2,0.5,1,2,5,10)) + # limits=c(ymin,ymax),
  scale_x_continuous(limits=c(0,1)) +
  guides(col=guide_legend(title="HR Endpoint 2")) + xlab('Correlation') + ylab('ARE') +
  theme(legend.position="bottom",
        axis.text=element_text(size=12),
        axis.title=element_text(size=15,face="bold"))  
ggsave('RPlot_ARET2E.png')  
