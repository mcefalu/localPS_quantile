library(quantreg)
library(splines)

lprq <- function(x, y, h, m=50 , tau=.5) {
   xx <- seq(min(x),max(x),length=m)
   fv<-xx
   dv<-xx
   for(i in 1:length(xx)) {
      z <- x - xx[i]
      wx <- dnorm(z/h)
      r <- rq(y~z, weights=wx, tau=tau, ci=FALSE)
      fv[i] <- r$coef[1.]
      dv[i] <- r$coef[2.]
   }
   list(xx = xx, fv = fv, dv = dv)
}

expit= function(x) exp(x)/(1+exp(x))
n = 1000
X = sort(rnorm(n))
ps = expit(-1+X)
A = rbinom(n,1,ps)
Y = rnorm(n,A+X)
ps = glm(A~X,family='binomial')$fitted

lm(Y~A)
lm(Y~X+A)


plot(ps, Y, xlab = "milliseconds", ylab = "acceleration", type = "n")
points(ps, Y, cex = 0.75)

XX <- model.matrix(Y ~ bs(ps, df = 15)) 


fit <- function(x) rq(Y ~ bs(ps, df = 15), tau = x)

accel.fit <- function(x) XX %*% fit(x)$coef

out=rowSums(accel.fit(seq(0,1,length=n)))/n

lines(ps, out , col='blue')


M=lprq(Y,ps,1)
lines(M$xx,M$fv)


M = 200
beta = matrix(0,M,3)
for (i in 1:M){
   X = sort(rnorm(n))
   ps = expit(-1+X)
   A = rbinom(n,1,ps)
   Y = rnorm(n,A+X-.5*A*X^2)
   ps = glm(A~X,family='binomial')$fitted
   w = A/ps + (1-A)/(1-ps)
   
   
   XX <- model.matrix(Y ~ bs(ps, df = 15)+A) 
   fit <- function(x) {
      K = dnorm( (
      rq(Y ~ bs(ps, df = 15)+A, tau = x)
   }
   fit2 <- lprq(ps,Y,1)
   lines(fit2$xx,fit2$fv)
   
   accel.fit <- function(x) fit(x)$coef["A"]
   
   out=sapply(seq(0,1,length=n),accel.fit)
   plot(ps, Y, xlab = "ps", ylab = "outcome", type = "n")
   points(ps, Y, cex = 0.75)
   lines(seq(0,1,length=n), out , col='blue')
   
   beta[i,1]=lm(Y~A+X)$coef["A"]
   beta[i,2]=mean(out)
   beta[i,3]=lm(Y~A,weights=w)$coef["A"]
}
beta = beta[1:135,]

   D = as.vector(beta)
   DD = rep(1:ncol(beta),each=nrow(beta))
   boxplot(D~DD)

colMeans(beta)
apply(beta,2,sd)

