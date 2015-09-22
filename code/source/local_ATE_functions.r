############################
###    LoWePS functions  ### 
##
## Recent modifications: Add necessary calculations so that the marginal
##                       quantile function is correctly calculated.
##
##                       Moved calculations for figures into my.localQ
##                       This allows us to not save the quantile regression fits.
## Date: 09/08/15
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

my.localQ <- function(y,a,ps,m=100,m.q=1000,h=0.1,boot=F,K.boot=100){
  if (is.numeric(m)){
    # sequence of PS from min to max PS with length m
    min.ps = max(min(ps[a==1])-h/2,0,min(ps[a==0]))
    max.ps = min(max(ps[a==0])+h/2,1,max(ps[a==1]))
    grid <- seq(min.ps,max.ps,length=m)
    index.ps <- 1:m
  }else{
    # use observed PS instead
    min.ps = max(min(ps[a==1])-h/2,0,min(ps[a==0]))
    max.ps = min(max(ps[a==0])+h/2,1,max(ps[a==1]))
    index.ps = which((ps<=max.ps) & (ps>=min.ps))
    grid <- ps[index.ps]
    m <- length(grid)
  }
  # grid of quantiles for collapsing results 
  q.grid = seq(0,1,length=m.q)
  n <- length(ps)
  n0 <- sum(1-a)
  n1 <- sum(a)
  
  # bandwidth -- currently specified by user
  #h = diff(range(grid))/h
  
  # some null vectors
  delta.mean = NULL
  q.heatmap = qY0.heatmap = qY1.heatmap = NULL
  delta.heatmap = NULL
  delta.q = matrix(NA,m.q,m)
  ps.heatmap = NULL
  K.ps <- K.ps1 <- K.ps0 <- numeric(m)
  for (i in 1:m){
      # weights
      K = dnorm( (ps-grid[i])/h )
      # exploratory: zero some weights out and standardize
      K[ K < dnorm(1) ] = 0 
      K <- K / (1-2*dnorm(-1))
      K.ps[i] = sum(K)/n/h/m
      K.ps0[i] = sum(K[a==0])/n0/h/m
      K.ps1[i] = sum(K[a==1])/n1/h/m
      
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
      
      ############################################
      ## add estimates of the marginal quantile ##
      # delta.q has the quantiles in the rows
      # !!!!!!!!!!!!!!! I DONT THINK THIS IS DOING WHAT I WANT
      delta.q[,i] = out$sol["a",findInterval(q.grid,out$sol["tau",])]
    
      #######################################
      # add estimates and PS to the output ##
      delta.heatmap = c(delta.heatmap,out$sol["a",])
      ps.heatmap = c(ps.heatmap,rep(grid[i],length(out$sol["a",])))
  }
  delta.mean.ub = delta.mean.lb = delta.mean.se = NULL
  if (boot){
    delta.mean.boot = matrix(NA,K.boot,m)
    message(paste("Performing",K.boot,"bootstraps:"))
    pb <- txtProgressBar(min=1,max=K.boot,style=3,char="=",width=50)
    for (i in 1:K.boot){
      setTxtProgressBar(pb, i)
      index = sample(n , n , replace=T)
      out = my.localQ(y=y[index],a=a[index],ps=ps[index],m=m,h=h)
      delta.mean.boot[i,] = out$delta.mean
    }
    delta.mean.lb = apply(delta.mean.boot,2,quantile,0.025)
    delta.mean.ub = apply(delta.mean.boot,2,quantile,0.975)
    delta.mean.se = apply(delta.mean.boot,2,sd)
  }
  # return the fits, the PS grid, the bandwidth, the outcomes , and the treatments
  return(list(ps=grid,index.ps=index.ps,q=q.grid,K.ps=K.ps,K.ps1=K.ps1,K.ps0=K.ps0,h=h,y=y,a=a,q.heatmap=q.heatmap,qY0.heatmap=qY0.heatmap,qY1.heatmap=qY1.heatmap,ps.heatmap=ps.heatmap,delta.heatmap=delta.heatmap,delta.mean=delta.mean,delta.mean.ub=delta.mean.ub,delta.mean.lb=delta.mean.lb,delta.mean.se=delta.mean.se,delta.q=delta.q))
}



###########################################
### plot mean as a function of PS based ###
### on local quantile regressions       ###
##
## fit is "out" from my.localQ
## add option adds plot to previous plot
## everything else is passed through to plot

plot.Qmean <- function(fit,add=F,boot=F,ps.dist=F,...){
  if (ps.dist){
    layout( rbind(matrix(2,4,5),1) )
    boxplot(fit$ps~fit$a)
  }
  if (!boot){
    if (add){
      index = order(fit$ps)
      lines(y=fit$delta.mean[index],x=fit$ps[index],...)
    }else{
      index = order(fit$ps)
      plot(y=fit$delta.mean[index],x=fit$ps[index],type='l',...)
    }
  }else{
    index = order(fit$ps)
    plot(fit$delta.mean[index],x=fit$ps[index],type='l',ylim=c(min(fit$delta.mean.lb),max(fit$delta.mean.ub)),...)
    lines(y=fit$delta.mean.lb[index],x=fit$ps[index],lty='dashed',...)
    lines(y=fit$delta.mean.ub[index],x=fit$ps[index],lty='dashed',...)
  }
}




#######################################
### plot marginal quantile function ###
##
## fit is "out" from my.localQ
## Ytype defines how we plot the x-axis.

plot.Qavg <- function(fit,Ytype="quantile",...){
  # i updated this to weight based on PS dist
  # but my implementation makes plotting the actual
  # Ys difficult
  
  #   if (Ytype=="observedY1"){
  #     q = fit$qY1.heatmap
  #   }
  #   if (Ytype=="observedY0"){
  #     q = fit$qY0.heatmap
  #   }
  #   if (Ytype=="quantile"){
  #     q = fit$q.heatmap
  #   }
  qs=fit$delta.q%*%fit$K.ps
  #tapply(fit$delta.heatmap,round(q,2),sum)
  plot(y=(qs),x=fit$q,...)
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



################################
### bootstrap for my.localQ ###
##
## y is outcome
## a is treatment
## ps is propensity score
## m is number of points on a grid for the PS

boot.localPSreg <- function(K=100,y,a,ps,m=100,h=0.1){
    delta.mean = matrix(NA,K,m)
    n = length(y)
    for (i in 1:K){
      index = sample(n , n , replace=T)
      out = my.localQ(y=y[index],a=a[index],ps=ps[index],m=m,h=h)
      delta.mean[i,] = out$delta.mean
    }
    ps.grid = out$ps
    return(list(delta=delta.mean,ps=ps.grid))
}
    