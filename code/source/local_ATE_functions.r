############################
###    LoWePS functions  ### 
##
## Last modification: add comments / clean code
## Date: 08/28/15
## By: M Cefalu

###############################
### load required packages  ### 
library(quantreg)
library(Hmisc)
library(KernSmooth)
library(survey)



####################################
### my local quantile regression ###
##
## y is outcome
## a is treatment
## ps is propensity score
## m is number of points on a grid for the PS

my.localQ <- function(y,a,ps,m=100,h=0.1){
  # sequence of PS from min to max PS with length m
  grid <- seq(min(ps),max(ps),length=m)
  # bandwidth -- currently specified by user
  #h = diff(range(grid))/h
  
  # out holds the quantile regression fits for each point in the PS grid
  out = list()
  for (i in 1:m){
    # weights
    K = dnorm( (ps-grid[i])/h )
    # exploratory: zero some weights out
    K[ K < dnorm(1) ] = 0 
    
    ######################################
    ## fit weighted quantile regression ##
    # tau is the quantile to be estimated
    # tau>1 estimates solution for all values of tau in (0,1)
    out[[i]] = rq(y~a,tau=3 , weights=K)
  }
  # return the fits, the PS grid, the bandwidth, the outcomes , and the treatments
  return(list(out=out,ps=grid,h=h,y=y,a=a))
}




###########################################
### plot mean as a function of PS based ###
### on local quantile regressions       ###
##
## out is "out" from my.localQ
## add option adds plot to previous plot
## everything else is passed through to plot

plot.Qmean <- function(fit,add=F,...){
  # delta saves estiamtes of the mean
  delta = NULL
  # loop over grid of PS used in my.localQ
  for (i in 1:length(fit$out)){
    # add estimates of the mean
    # the mean is just the integral of the quantile function
    # for each PS, the quantile function is a step function
    # so we just add up the estimates, weighting by the width between quantile points
    delta = c(delta,sum(diff(fit$out[[i]]$sol["tau",])*fit$out[[i]]$sol["a",-ncol(fit$out[[i]]$sol)]))
  }
  # now plot
  if (add){
    lines(y=delta,x=fit$ps,...)
  }else{
    plot(y=delta,x=fit$ps,type='l',...)
  }
}




#######################################
### plot marginal quantile function ###
##
## fit is "out" from my.localQ
## Ytype defines how we plot the x-axis.

plot.Qavg <- function(fit,Ytype="observedY0",...){
  # set outcome and treatment
  Y=fit$y
  X=fit$a
  # delta will hold effect estimates
  delta = NULL
  # q are the quantiles of Y(0) or Y(1), or just the quantile itself
  q = NULL
  # x is the propensity score
  x = NULL
  N = length(Y)
  
  # loop over the grid of propensity scores used in my.localQ
  for (i in 1:length(fit$out)){
    # weights used for the i-th quantile reg
    K = fit$out[[i]]$weights
    # check ytype -- it changes what we save in q
    if (Ytype=="observedY1"){
      temp = wtd.quantile(Y[X==1],probs=fit$out[[i]]$sol["tau",],weights=K[X==1])
      q = c(q,temp)
    }
    if (Ytype=="observedY0"){
      temp = wtd.quantile(Y[X==0],probs=fit$out[[i]]$sol["tau",],weights=K[X==0])
      q = c(q,temp)
    }
    if (Ytype=="observedY"){
      temp = wtd.quantile(Y,probs=fit$out[[i]]$sol["tau",],weights=K)
      q = c(q,temp)
    }
    if (Ytype=="quantile"){
      q = c(q,fit$out[[i]]$sol["tau",])
    }
    # add estimates and PS to the output
    delta = c(delta,fit$out[[i]]$sol["a",])
    x = c(x,rep(fit$ps[i],length(fit$out[[i]]$sol["a",])))
  }
  # might need to weight this based on observed PS dist?
  qs=tapply(delta,round(q,2),mean)
  plot(x=rownames(qs),y=qs,...)
}





###############################
###  heat map for my.localQ ###
##
## fit is "out" from my.local Q
## Ytype defines how we plot the y-axis.

heatmap.localQ <- function(fit,Ytype="observedY0",...){
  # set outcome and treatment
  Y=fit$y
  X=fit$a
  # check that ytype is valid
  if ( !(Ytype%in%c("observedY1","observedY0","observedY","quantile")) ){
    print("Invalid Ytype")
    return(NULL)
  }
  # effect estimates
  delta = NULL
  # quantiles of Y0
  q = NULL
  # propensity scores
  x = NULL
  N = length(Y)
  # this code is the same as from plot.Qavg
  for (i in 1:length(fit$out)){
    # weights used for the i-th quantile reg
    K = fit$out[[i]]$weights
    # check ytype -- it changes what we save in q
    if (Ytype=="observedY1"){
      temp = wtd.quantile(Y[X==1],probs=fit$out[[i]]$sol["tau",],weights=K[X==1])
      q = c(q,temp)
      ylab="Y1"
    }
    if (Ytype=="observedY0"){
      temp = wtd.quantile(Y[X==0],probs=fit$out[[i]]$sol["tau",],weights=K[X==0])
      q = c(q,temp)
      ylab = "Y0"
    }
    if (Ytype=="observedY"){
      temp = wtd.quantile(Y,probs=fit$out[[i]]$sol["tau",],weights=K)
      q = c(q,temp)
      ylab = "Y"
    }
    if (Ytype=="quantile"){
      q = c(q,fit$out[[i]]$sol["tau",])
      ylab = "quantile of Y(0)"
    }
    # add estimates and PS to the output
    delta = c(delta,fit$out[[i]]$sol["a",])
    x = c(x,rep(fit$ps[i],length(fit$out[[i]]$sol["a",])))
  }
  # plot! 
  index = q<Inf
  col = heat.clr(delta[index])
  plot(x=x[index],y=q[index],col=col,ylab=ylab,...)
}





##########################
### colors for heatmap ###
##
## grid is the values to be plotted
## color is 3 color system. First color is negative, second zero, third positive

heat.clr = function(grid,color=c('red','white',"blue"),...){
  clr1 = as.vector(grid)#-mean(grid)
  clr.final = numeric(length(clr1))
  c.max = max(abs(clr1))
  c.min = min(abs(clr1))
  clr2 = (clr1-c.min)/(c.max-c.min)
  # for positives
  index = clr1>0
  shade = colorRamp(color[2:3])(clr2[index])
  color.hex = rgb(shade,max=255)
  clr.final[index] = color.hex
  # for negatives
  index = clr1<=0
  shade = colorRamp(color[2:1])(-clr2[index])
  color.hex = rgb(shade,max=255)
  clr.final[index] = color.hex
  return(clr.final)
}




###########################
### local ps regression ###
## 
## Y is outcome
## X is treatment
## ps is propensity score
## h is bandwidth

localPSreg <- function(Y,ps,X,h=.1){
  V=beta=numeric(length(Y))
  CI = matrix(0,length(Y),2)
  # loop over observations and do a weighted regression around that obs
  for (i in 1:length(Y)){
    # bandwidth -- user specified for now
    #h = diff(range(ps))/15
    # weights
    K = dnorm( (ps[i]-ps)/h )/dnorm(0) 
    # fit the model
    dta = data.frame(Y=Y,X=X,K=K,ps=ps)
    D <- svydesign(id = ~1, weights = ~K, data=dta)
    M = svyglm(Y~X,D)
    # save output
    beta[i] = M$coef[2]
    V[i] = vcov(M)[2,2]
    CI[i,] = beta[i] + c(-1.96,1.96)*sqrt(V[i])
  }
  return(list(beta=beta,CI=CI,V=V))
}



#####################################
### plot mean based on localPSreg ###
##
## fit is output of localPSreg
## CI specifies whether or not to plot CIs
## add specifies whether or not to add to existing figure
## everything else passed to plot

plot.localPS <- function(out,CI=F,add=F,...){
  # this code is broken!!!
  index = order(dta$ps)
  ps = dta$ps[index]
  beta.localPS <- out.localPS$beta[index]
  if (CI){
    CI.localPS <- out.localPS$CI[index,]
    ylim = range(c(beta.localPS,CI.localPS))
    if (add){
      lines(ps,beta.localPS,ylim=ylim,...)
      lines(ps,CI.localPS[,1],lty="dashed",col="grey")
      lines(ps,CI.localPS[,2],lty="dashed",col="grey")
    }else{
      plot(ps,beta.localPS,type='l',ylim=ylim,...)
      lines(ps,CI.localPS[,1],lty="dashed",col="grey")
      lines(ps,CI.localPS[,2],lty="dashed",col="grey")
    }
  }else{
    if (add){
      lines(ps,beta.localPS,...)
    }else{
      plot(ps,beta.localPS,type='l',xlim=c(0,1),...)
    }
  }
}
