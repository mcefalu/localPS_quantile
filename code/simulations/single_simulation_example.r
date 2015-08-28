###########################################
### code that simulates single dataset  ### 
##
## Last modification: add comments / clean code
## Date: 08/28/15
## By: M Cefalu

rm(list=ls())
source("code/source/local_ATE_functions.r")
source("code/source/simulation_functions.r")

# generate data
set.seed(14322)
dta = gen.dta(n=1000,a=c(-1,1),b=c(0,1),bx.mult=-1,bx.add=1)

# fit the PS model
m=glm(X~C,family='binomial',data=dta)
dta$ps=m$fitted

# regression adjustment for PS
lm(Y~X+ps,data=dta)

# local PS reg and local Q reg
out <- localPSreg(Y=dta$Y,X=dta$X,ps=dta$ps)
out <- my.localQ(y=dta$Y,a=dta$X,ps=dta$ps)

# plot estimates of mean, conditional on PS
# the plot.localPS function is currently broken!
#plot.localPS(fit=out,CI=T,add=F)
plot.Qmean(fit=out,add=F,col='blue')
  
# plot truth -- the f() function is no longer correct.
#lines(sort(ps),f(sort(ps),a=a,b=b,bx=bx),col='red')

# my heat map 
heatmap.localQ(fit=out,Ytype="quantile",pch=16,xlab="PS")

## plot effect at each quantiles -- marginalized across PS
plot.Qavg(fit=out,Ytype="observedY0",type='l',col='blue',xlab="Y0",ylab="Effect")

