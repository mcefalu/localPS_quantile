############################
### Simulation functions ### 
##
## Last modification: new file
## Date: 08/28/15
## By: M Cefalu


# expit function
expit <- function(x) exp(x)/(1+exp(x))

# generate data
gen.dta <- function(n,a=c(0,0),b=c(0,0),bx.mult=0,bx.add=0,s2=1){
  # generate covariate
  C = rnorm(n)
  # generate treatment
  X = rbinom(n,1,expit(cbind(1,C)%*%a))
  # coefs for outcome
  Y0 = rnorm(n,cbind(1,C)%*%b,sqrt(s2))
  Y1 = Y0*exp(bx.mult) + bx.add
  # observed Y
  Y = X*Y1 + (1-X)*Y0
  return(data.frame(Y=Y,X=X,C=C,Y0=Y0,Y1=Y1))
}

# transform from PS to treatment effect
# NOT UP TO DATE!!!
f = function(p,a,b,bx) {
  # x is really C
  x = (log(p/(1-p))-a[1])/a[2]
  Y1 = g(cbind(1,x,1,x^2)%*%b)*exp(bx)
  Y0 = g(cbind(1,x,0,0)%*%b)
  return(Y1-Y0)
}

# function to calculate ATE
# NOT UP TO DATE!!!!
ate = function(b) {
  n = 1000000
  C = rnorm(n)
  Y1 = g(cbind(1,C,1,C)%*%b)
  Y0 = g(cbind(1,C,0,0)%*%b)
  return(mean(Y1-Y0))
}
