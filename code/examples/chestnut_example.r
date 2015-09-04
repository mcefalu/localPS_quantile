rm(list=ls())
source("code/source/local_ATE_functions.r")
library(twang)

# read data from BA
dta=read.csv("data/raw/data_from_BA.csv")

# fit PS model using GBM
ps1 <- ps(atm ~ age + female + female + race4g + sfs + sps + sds + ias + ces + eps + imds + bcs + prmhtx,
          data = dta,
          n.trees = 5000,
          interaction.depth = 3,
          shrinkage = 0.01,
          perm.test.iters = 0,
          stop.method = c("es.max"),
          estimand = "ATE",
          sampw = NULL,
          verbose = FALSE
)

# check balance
bal.table(ps1)

# Also – here is a list of outcomes to consider
# sfs8p12 = SFS at 12-months
# spsm12 = SPS at 12-months
# eps7p12 = EPS at 12-months
# sdsm12 = SDS at 12-months
# Maxce12 = (indicator of being institutionalized in last 90 days

# save the ps
dta$ps = ps1$ps[,1]

# fit LoWePS-QR
# h=dpill(x=lalonde$ps, y=lalonde$re78)
out.sfs <- my.localQ(y=dta$sfs8p12,a=dta$atm,ps=dta$ps,h=0.1)
out.sps <- my.localQ(y=dta$spsm12,a=dta$atm,ps=dta$ps,h=0.1)
out.eps <- my.localQ(y=dta$eps7p12,a=dta$atm,ps=dta$ps,h=0.1)
out.sds <- my.localQ(y=dta$sdsm12,a=dta$atm,ps=dta$ps,h=0.1)

save.image("data/results/chesnut_results.Rdata")

# generate plots!
#heatmap.localQ(fit=out.sfs , Ytype="quantile" , pch=16 , xlab="PS")
#plot.Qmean(fit=out.sfs , add=F , col='blue' , xlab="PS" , ylab="Effect")
#plot.Qavg(fit=out.sfs , Ytype="quantile" , type='l' , col='blue' , xlab="quantile of Y(0)" , ylab="Effect")











