install.packages("forecast")
library(forecast)

install.packages("tseries")
library(tseries)

install.packages("urca")
library(urca)

seed <- 423
set.seed(seed)

frequency <- 5
lag_max_acf <- 400
lag_max_pacf <- 400
diff <- 1

xau_df <- read.csv("XAU.csv", header = TRUE, stringsAsFactors = FALSE)
close_all <- xau_df$Close.Last
close_420 <- tail(close_all, 420)
y_test  <- close_420[1:20]
y_train <- close_420[21:420]


ts_train <- ts(y_train,frequency = frequency)


length(y_train)    
length(y_test)     
plot(ts_train, main = "Training series (last 400 values of Close.Last)")

cfs<-function(ts_train,l1,l2){
  acf(ts_train,
      main = "ACF",
      lag.max=l1)          
  

  pacf(ts_train,
       main = "PACF",
       lag.max=l2)         
}

cfs(ts_train,lag_max_acf,lag_max_pacf)

stationarity_check <- function(ts_train){
  adf_result <- adf.test(ts_train,       
                         alternative = "stationary",
  )
  print(adf_result)
  

  kpss_level <- kpss.test(ts_train, null = "Level")
  print(kpss_level)

  kpss_trend <- kpss.test(ts_train, null = "Trend")
  print(kpss_trend)
}

stationarity_check(ts_train)


lambda <- BoxCox.lambda(ts_train, lower = -3, upper = 3)
cat("Estimated lambda:", lambda, "\n")


ts_train_bc <- BoxCox(ts_train, lambda)
plot(ts_train_bc, main = paste0("Box-Cox (λ=", round(lambda,2), ")"))

stationarity_check(ts_train_bc)



ts_d1 <-diff(ts_train,differences=diff)

plot(ts_d1,
     main  = expression(paste("First Difference of ", ts[train])),
     ylab  = expression(Delta * x[t]),
     xlab  = "Time",
     type  = "l",      
     lwd   = 2,        
     col   = "blue")    


abline(h = 0, lty = 2, col = "gray")
cfs(ts_d1,lag_max_acf,lag_max_pacf)

stationarity_check(ts_d1)

max_order <- 4
orders <- list()
orders <- lapply(1:max_order, function(i) c(i, 0, i))

#https://www.w3schools.com/r/r_data_frames.asp
results <- data.frame(
  model = character(0),
  AIC   = numeric(0),
  BIC   = numeric(0),
  SS    = numeric(0),
  stringsAsFactors = FALSE
)

for (o in orders) {
  fit <- tryCatch(
    arima(ts_d1, order = o),
    error = function(e) NULL
  )
  if (!is.null(fit)) {
    results <- rbind(
      results,
      data.frame(
        model = paste0("(", paste(o, collapse = ","), ")"),
        AIC   = round(fit$aic, 3),
        BIC   = round(BIC(fit), 3),
        SS    = round(mean(fit$residuals^2), 3),
        stringsAsFactors = FALSE
      )
    )
  }
}

print(results)

chosen_order <- 1
orders <- list()
for(i in 0:chosen_order) {
  orders[[  2*(i+1) -1  ]] <- c(i, 0, chosen_order)
  orders[[  2*(i+1)   ]] <- c(chosen_order, 0, i)
  
}


#https://www.w3schools.com/r/r_data_frames.asp
results <- data.frame(
  model = character(0),
  AIC   = numeric(0),
  BIC   = numeric(0),
  SS    = numeric(0),
  stringsAsFactors = FALSE
)


for (o in orders) {
  fit <- tryCatch(
    arima(ts_d1, order = o),
    error = function(e) NULL
  )
  if (!is.null(fit)) {
    results <- rbind(
      results,
      data.frame(
        model = paste0("(", paste(o, collapse = ","), ")"),
        AIC   = round(fit$aic, 3),
        BIC   = round(BIC(fit), 3),
        SS    = round(mean(fit$residuals^2), 3),
        stringsAsFactors = FALSE
      )
    )
  }
}

print(results)


model_0_1_1 <-arima(ts_train, order=c(0,1,1))
tsdiag(model_0_1_1)


model_0_1_1$coef
model_0_1_1$sigma2


rs <- model_0_1_1$residuals
stdres <- rs/sqrt(model_0_1_1$sigma2)
qq <- qqnorm(stdres, plot.it = FALSE)


plot(qq$x, qq$y,
     main = "QQ Plot of Standardized Residuals\nwith Regression Line",
     xlab = "Theoretical Quantiles",
     ylab = "Sample Quantiles",
     pch  = 20)


fit <- lm(qq$y ~ qq$x)
abline(fit, col = "magenta", lwd = 1)



abline(0, 1, col = "blue", lty = 2)

legend("topleft",
       legend = c("Sample Quantiles", "Regression Line", "First Diagonal Line"),
       col    = c("black",            "magenta",             "blue"),
       pch    = c(20,                 NA,                NA),
       lty    = c(NA,                 1,                 2),
       lwd    = c(NA,                 2,                 1),
       bg     = "white")


auto_fit <- auto.arima(
  ts_train,
  max.p = 5,
  max.q = 5,
  max.P = 3,
  max.Q = 3,
  max.d = 2,
  max.D = 2,
  max.order = 20,
  
  stepwise   = FALSE,   
  approximation = FALSE,
  trace      = TRUE     
)


f011  <- predict(model_0_1_1, n.ahead = 20)
fauto <- predict(auto_fit,      n.ahead = 20)

pred011 <- as.numeric(f011$pred)
predauto <- as.numeric(fauto$pred)


actual <- c(y_train[1:20], y_test)   

plot(actual,
     type = "l",
     lwd  = 2,
     col  = "black",
     xlab = "Time (predictions start at 21)",
     ylab = "Value",
     main = "20-Step Forecasts vs. Actuals")


lines(x = 21:40, y = pred011,
      col = "red",  lwd = 2, lty = 1)


lines(x = 21:40, y = predauto,
      col = "blue", lwd = 2, lty = 2)

legend("topleft",
       legend = c("Actual (hold-out + test)",
                  "ARIMA(0,1,1) Forecast",
                  "auto.arima Forecast"),
       col    = c("black", "red",    "blue"),
       lty    = c(1,       1,        2),
       lwd    = 2,
       bty    = "n")







