############################
###    LoWePS functions  ### 
##
## Last modification: moved calculations for figures into my.localQ
##                    This allows us to not save the quantile regression fits.
## Date: 09/04/15
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
  # out = list()
  delta.mean = NULL
  q.heatmap = qY0.heatmap = qY1.heatmap = NULL
  delta.heatmap = NULL
  ps.heatmap = NULL
  for (i in 1:m){
    # weights
    K = dnorm( (ps-grid[i])/h )
    # exploratory: zero some weights out
    K[ K < dnorm(1) ] = 0 
    
    ######################################
    ## fit weighted quantile regression ##
    # tau is the quantile to be estimated
    # tau>1 estimates solution for all values of tau in (0,1)
    out = rq(y~a,tau=3 , weights=K)
    
    ###############################
    ## add estimates of the mean ##
    # the mean is just the integral of the quantile function
    # for each PS, the quantile function is a step function
    # so we just add up the estimates, weighting by the width between quantile points
    delta.mean = c(delta.mean,sum(diff(out$sol["tau",])*out$sol["a",-ncol(out$sol)]))
    
    ##########################################
    ## extract quantiles for use in heatmap ##
    q.heatmap = c(q.heatmap,out$sol["tau",])
    temp = wtd.quantile(y[a==1],probs=out$sol["tau",],weights=K[a==1])
    qY1.heatmap = c(qY1.heatmap,temp)
    temp = wtd.quantile(y[a==0],probs=out$sol["tau",],weights=K[a==0])
    qY0.heatmap = c(qY0.heatmap,temp)
  
    #######################################
    # add estimates and PS to the output ##
    delta.heatmap = c(delta.heatmap,out$sol["a",])
    ps.heatmap = c(ps.heatmap,rep(grid[i],length(out$sol["a",])))
  }
  # return the fits, the PS grid, the bandwidth, the outcomes , and the treatments
  return(list(ps=grid,h=h,y=y,a=a,q.heatmap=q.heatmap,qY0.heatmap=qY0.heatmap,qY1.heatmap=qY1.heatmap,ps.heatmap=ps.heatmap,delta.heatmap=delta.heatmap,delta.mean=delta.mean))
}




###########################################
### plot mean as a function of PS based ###
### on local quantile regressions       ###
##
## fit is "out" from my.localQ
## add option adds plot to previous plot
## everything else is passed through to plot

plot.Qmean <- function(fit,add=F,...){
  if (add){
    lines(y=fit$delta.mean,x=fit$ps,...)
  }else{
    plot(y=fit$delta.mean,x=fit$ps,type='l',...)
  }
}




#######################################
### plot marginal quantile function ###
##
## fit is "out" from my.localQ
## Ytype defines how we plot the x-axis.

plot.Qavg <- function(fit,Ytype="quantile",...){
  # might need to weight this based on observed PS dist?
  if (Ytype=="observedY1"){
    q = fit$qY1.heatmap
  }
  if (Ytype=="observedY0"){
    q = fit$qY0.heatmap
  }
  if (Ytype=="quantile"){
    q = fit$q.heatmap
  }
  qs=tapply(fit$delta.heatmap,round(q,2),mean)
  plot(x=rownames(qs),y=qs,...)
}





###############################
###  heat map for my.localQ ###
##
## fit is "out" from my.local Q
## Ytype defines how we plot the y-axis.

heatmap.localQ <- function(fit,Ytype="quantile",...){
  if (Ytype=="observedY1"){
    q = fit$qY1.heatmap
  }
  if (Ytype=="observedY0"){
    q = fit$qY0.heatmap
  }
  if (Ytype=="quantile"){
    q = fit$q.heatmap
  }
  index = q<Inf
  col = heat.clr(fit$delta.heatmap[index])
  plot(x=fit$ps.heatmap[index],y=q[index],col=col,...)
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
