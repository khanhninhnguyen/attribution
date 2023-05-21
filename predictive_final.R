# Final predictive rule prog

# output: last predtive rule and prediction results 
# input: test result, name of series used for training: 
# G-E, G-G', G-E',E-E', G'-E', G'-E, remove them by name of remove variable 
# prob for each configuration, significance.level, name.version


# load package and paths --------------------------------------------------

if (!require(devtools)) install.packages("devtools")
devtools::install_github("yanlinlin82/ggvenn")
library(ggrepel)
library(mclust)
library(forecast)
library(lmtest)
library(ggvenn)
library(fields)
library(MASS)
library(reshape2)
library(factoextra)
library(mvtnorm)
library(caret)
library(randomForest)
library(rpart)
library(rpart.plot)
library(class)
library(rsample)
library(tidyverse)

path_restest <- paste0(path_results,"attribution/predictive_rule/")
file_path_Results=paste0(path_results,'attribution/predictive_rule/')


# make the truth table and probability ------------------------------------

G=c(rep(1,9), rep(0,9),rep(-1,9),rep(0,9),rep(1,9),rep(-1,9))
E=c(rep(0,9),rep(-1,9),rep(0,9),rep(1,9),rep(-1,9),rep(1,9))
a=c(rep(0,3),rep(1,3),rep(-1,3))
Gp=rep(a,6)
Ep=rep(c(0,1,-1),18)
Y=data.frame(G=G,E=E,Gp=Gp,Ep=Ep)

Z=data.frame(GE=Y$G-Y$E,GGp=Y$G-Y$Gp,GEp=Y$G-Y$Ep,EEp=Y$E-Y$Ep,GpEp=Y$Gp-Y$Ep,GpE=Y$Gp-Y$E)
knitr::kable(head(Z))
List.names.tot <- colnames(Z)

Z.trunc <- Z 
Z.trunc[Z.trunc==-2]=-1
Z.trunc[Z.trunc==2]=1
knitr::kable(head(Z.trunc))

prob <- c(0.18225,0.010125,0.010125,0.010125,0.0005625,0.0005625,0.010125,0.0005625,0.0005625,
          0.010125,0.010125,0.18225,0.0005625,0.0005625,0.010125,0.0005625,0.0005625,0.010125,
          0.18225,0.010125,0.010125,0.010125,0.0005625,0.0005625,0.010125,0.0005625,0.0005625,
          0.010125,0.18225,0.010125,0.0005625,0.010125,0.0005625,0.0005625,0.010125,0.0005625,
          0.00225,0.00225,0.0405,0.010125,0.010125,0.00225,0.010125,0.010125,0.00225,0.00225,
          0.0405,0.00225,0.010125,0.00225,0.010125,0.010125,0.00225,0.010125)

keep.config <- c(1:3,6:15,17,19:24,26,28:30,33:40,43,46:49,52)
Y <- Y[keep.config,]
Z <- Z[keep.config,]
Z.trunc <- Z.trunc[keep.config,]
num.conf <- 1:length(keep.config)
row.names(Z.trunc) <- num.conf

prob <- prob[keep.config]
prob <- prob/sum(prob)

# remove the duplicated due to the selection of variable ------------------

rm.ind = which(list.name.test == remove.var)
sum(duplicated2(Z.trunc[,-rm.ind]))/2
rg.duplicate <- which(duplicated2(Z.trunc[,-rm.ind ])) 









