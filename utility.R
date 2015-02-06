## Richard D. Morey
## January 2015

## Included for some functions related to 
## arithmetic with logarithms
library(BayesFactor)

## Find tau such that the error rates are equal
calibrated_tau = Vectorize(function(d1, N, k = 10, NRshift = .1, tol = .01){
  # Find critical delta value for the error rates
  deltahat = crit_delta(d1,N)
  # Convert to critical t
  t = deltahat * sqrt(N)
  # Find the optimal tau, then move to the right a small amount
  # to start Newton-Raphson, because we want the non-trivial solution
  # if it exists (that is, tau != 0)
  starttau = optimal_tau(N,t) + NRshift
  # Convert tau to qz = 1 + N * tau^2 
  # for optimization 
  startqz = 1 + N*starttau^2
  zero_q = newtonSteps(startqz, f = f2, fp = df2dq, k=k, N=N, t=t)
  if(is.nan(zero_q)){
    warning("Algorithm failed to converge (extreme values?)")
    return(NA)
  }
  # Convert back to tau
  tau = sqrt((zero_q - 1)/N)
  # check to make sure algorithm converged within tolerance
  lbf = logBF(t, N, tau)
  if( ( lbf > log(1+tol) ) | ( lbf < log(1 - tol) ) ){
    warning("Algorithm failed to converge to within ", tol*100, "% tolerance. (final BF was off by ",exp(lbf)*100,"%)")
  }
  return(tau)
},c("d1","N"))


## Calibration restriction holds when probDiff = 0
## Uses d_crit = t_crit/sqrt(N)
probDiff = Vectorize(function(d_crit, d1, N){
  logpnull = log(2) + pt(-d_crit*sqrt(N), N - 1, log=TRUE,lower.tail = TRUE)
  logpnull = BayesFactor:::logExpXminusExpY(0,logpnull)
  logpalt0 = pt(-d_crit*sqrt(N), N-1, d1*sqrt(N),log=TRUE) 
  logpalt1 = pt(d_crit*sqrt(N), N-1, d1*sqrt(N), log = TRUE, lower.tail=FALSE)
  logpalt = logpalt1 + log1p(exp(logpalt0-logpalt1))
  (logpalt - logpnull)^2
}, "d_crit")

## Finds the estimate of dhat that would
## lead to equal errors
crit_delta = Vectorize(function(d1,N){
  optimize( probDiff, c(0, d1), N = N, d1 = d1)$minimum
},c("d1","N"))

## Find the calibrated d1 such that
## the criterion has a given t statistic, for 
## a specific N
find_d1 = Vectorize(function(t,N){
  alpha = 2 * ( 1 - pt( t, N-1) )
  optimize(function(d1, N, t, alpha){
    (pt(t, N - 1, ncp = d1*sqrt(N)) - alpha)^2  
  }, interval = c(t/sqrt(N),2*t/sqrt(N)), N = N, t = t, alpha = alpha)$minimum
}, c("t","N"))

## (log) scaled information Bayes factor
logBF = Vectorize(function(t, N, tau){
  -N/2 * log(1 + t^2/(N-1)) + .5 * log(1 + N*tau^2) + N/2 * log(1 + t^2 / ( (1 + N * tau^2) * ( N - 1 ) ) )  
}, c("t","N","tau"))

## derivative of the (log) scaled information
## Bayes factor
dlogBFdt = Vectorize(function(t, N, tau){
  -N / (1 + t^2/(N-1)) * t/(N-1) + N / (1 + t^2 / ( (1 + N * tau^2) * ( N - 1 ) ) ) * t / ( (1 + N * tau^2) * ( N - 1 ) )
},"t")

## Utility fuction for Newton-Raphson
## solver, given generic functions
newtonSteps = function(x, f, fp, k=3, ...)
  ifelse(k==0, 
         x, 
         newtonSteps(x - f(x, ...)/fp(x, ...), f=f, fp=fp, k=k-1, ...)
  )  

## For a given qz = 1 + N*tau^2, N, and t, report the
## (log) scaled information Bayes factor
f2 = Vectorize(function(qz, N, t){
  tau = sqrt((qz - 1)/N)
  logBF(t, N, tau)
},"qz")

## Derivative of the (log) scaled information
## Bayes factor, as a function of qz = 1 + N*tau^2
df2dq = Vectorize(function(qz, N, t){
  .5/qz * ( 1 - (N/(N-1)) * t^2 / (qz + t^2/(N-1))) 
},"qz")

## Find tau such that the (log) Bayes factor is minimized
optimal_tau = function(N, t){
  ifelse(t>1,sqrt((t^2 - 1)/N),0)
}


################
## For plots
################


## For reproducing Figures in article;
## Find the critical t for a given (log) Bayes factor to
## be equal to 0, and report errors under null and "alt"
probs = function(tau, d1, N, k = 4){
  t_crit = newtonSteps(1, f = logBF, fp = dlogBFdt, k=k, N=N, tau=tau)
  
  # Probability of B<1 (|t|<|t_c|) under null
  pnull = 1 - 2*pt(-t_crit, N - 1)
  # Probability of B>1 (|t|>|t_c|) under alternative
  palt = pt(-t_crit, N-1, d1*sqrt(N)) + (1-pt(t_crit, N-1, d1*sqrt(N)))
  c(t_crit,pnull,palt)
}


## Create distribution plot to show equal errors
distPlot = function(d1, N, min_d = -1, max_d=2, extrad = NULL){
  cd = crit_delta(d1, N)
  
  dd = seq(min_d, max_d, len=400)
  
  nully = dt(dd*sqrt(N),N-1)/sqrt(N)
  alty =  dt(dd*sqrt(N),N-1,ncp=d1 * sqrt(N))/sqrt(N)   
  
  par(lwd=2)
  plot(dd, nully, typ='l', col="red", ylab="", xlab= "Observed effect size", yaxt='n')
  lines(dd, alty, col="blue",lty=2)
  abline(v = c(-1,1)*cd, lty=4, col="gray")
  abline(h=0, col="gray")
  
  axis(3, at = c(-1,1)*cd, lab=c("BF=1","BF=1"))
  #mtext("Bayes factor", 3, 2.5, adj=.5, cex=1.2)
  
  mtext("Density", 2, 1, padj=.5, las=0, cex = 1.2)
  if(!is.null(extrad)){
    alty2 =  dt(dd*sqrt(N),N-1,ncp=extrad * sqrt(N))/sqrt(N)      
    lines(dd, alty2, col="blue",lty=3)
  }
}


