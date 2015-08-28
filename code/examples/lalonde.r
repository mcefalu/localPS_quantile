rm(list=ls())
source("code/source/local_ATE_functions.r")
library(twang)
data(lalonde)

# fit the PS model
ps.lalonde <- ps(treat ~ age + educ + black + hispan + nodegree +
married + re74 + re75,
data = lalonde, stop.method=c("ks.max"),
estimand = "ATE",verbose=FALSE)

# save the ps
lalonde$ps = ps.lalonde$ps[,1]

# fit LoWePS-QR
h=dpill(x=lalonde$ps, y=lalonde$re78)
out <- my.localQ(y=lalonde$re78,a=lalonde$treat,ps=lalonde$ps,h=0.1)

# generate plots!
heatmap.localQ(fit=out , Ytype="quantile" , pch=16 , xlab="PS")
plot.Qmean(fit=out , add=F , col='blue' , xlab="PS" , ylab="Effect")
plot.Qavg(fit=out , Ytype="quantile" , type='l' , col='blue' , xlab="quantile of Y(0)" , ylab="Effect")









