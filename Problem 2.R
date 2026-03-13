
install.packages("tseries")
library(tseries)

n_obs <- 200
seed=420
alpha <- 0.05
lag_max_acf <- 200
lag_max_pacf <- 200
p_order <- 2
ar <- c(1,-0.9)
lag_max_emp_pacf <- 200

set.seed(seed)
x<-arima.sim(model = list(order=c(p_order,0,0), ar=ar),
             n=n_obs)

plot.ts(x,main="My time series")
acf_res <-acf(x, lag.max=lag_max_acf, na.action = na.pass)
pacf(x, lag.max=lag_max_pacf, na.action = na.pass)
threshold <- qnorm(1-alpha/2) / sqrt(n_obs)


x_bar <- mean(x)
dev <- x - x_bar

cov <-function(h){
  return(sum ( dev[(h+1):n_obs]*dev[1:(n_obs-h)] ) /n_obs  )
}
gamma <- sapply(0:(n_obs-1),cov)

gamma_0=gamma[1]
gamma <- gamma[2:n_obs]

lag_vec <- as.numeric(acf_res$acf)

plot(c(gamma_0,gamma)/gamma_0,
     type = "h",                   
     lwd  = 1,                     
     xlab = "Lag",
     ylab = "Empirical Autocorrelation",
     main = "My (Empirical) ACF")
abline(h = threshold, col = "magenta", lty = 2)
abline(h = -threshold, col = "magenta", lty = 2)
abline(h = 0, col = "black")

plot(lag_vec-c(gamma_0,gamma)/gamma_0,
     type = "h",                   
     lwd  = 1,                    
     xlab = "Lag",
     ylab = "Diff",
     main = "R ACF - My ACF")

abline(h = 0, col = "black")


phi_nn<-gamma[1]/gamma_0
phi_vec<-phi_nn
V<-(gamma_0**2 - gamma[1]**2)/gamma_0

pacf_emp <- phi_nn
for(k in 2:(n_obs-1)){
  phi_nn <- (gamma[k] - sum( phi_vec*gamma[(k-1):1] )  ) / V
  phi_vec <- phi_vec - phi_nn*phi_vec[(k-1):1]
  phi_vec <- c(phi_vec,phi_nn)
  
  V <- V * (1-phi_nn**2)
  
  pacf_emp<-c(pacf_emp, phi_nn)
  
}

pacf(x, lag.max=50, na.action = na.pass)

plot(pacf_emp[1:lag_max_emp_pacf],
     type = "h",                   
     lwd  = 1,                     
     xlab = "Lag",
     ylab = "PACF",
     main = "My PACF")
abline(h = 0, col = "black")
abline(h = threshold, col = "magenta", lty = 2)
abline(h = -threshold, col = "magenta", lty = 2)

idxs <- which(abs(pacf_emp) >= threshold)
p   <- if (length(idxs)>0) max(idxs) else 0
p <- 1
for(i in length(idxs):1){
  if(i==1){
    p <- 1
    break
  }
  if(i==2){
    if(idxs[i]-1== idxs[i-1])
    {p <- 2
    break}
  }
  
  if(idxs[i]-1== idxs[i-1] && idxs[i-1]-1 == idxs[i-2])
  {p <- idxs[i]
  break}
  
}

print(p)