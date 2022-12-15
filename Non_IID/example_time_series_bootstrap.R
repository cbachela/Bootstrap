    
    
  # Code is taken from:
  # https://lbelzile.github.io/timeseRies/boostrap-methods-for-time-series.html

    
  
  #Counts are overdispersed
  library(boot)
  library(forecast)
  #Variance stabilizing transform
  sun <- 2*(sqrt(sunspot.year+1)-1)
  #QQ plot of the variance-transformed observations
  qqnorm(scale(sun), pty="s"); abline(a = 0, b = 1)
  
  
  #Alternative would be a Box-Cox transform
  # sun_bc <- BoxCox(window(sunspot.year, 1818), 
  #                     forecast::BoxCox.lambda(window(sunspot.year, 1818), method = "loglik"))
  # qqnorm(scale(sun_bc)); abline(a = 0, b = 1)
  # plot(sun_bc, main = "Average number of sunspots\n(Box-Cox transformed)", ylab = "")
  #apparent cycle of 11
  
  #Fit a time series model to transformed sunspot
  sun_ar_auto <- forecast::auto.arima(window(sun, 1930, 1979), max.q = 0, max.Q = 0, 
                                      max.d = 0, allowdrift = FALSE, max.D = 0, max.P = 0, 
                                      max.p = 25, max.order = 25, stepwise = FALSE, ic = "aic")
  p <- sun_ar_auto$arma[1]
  res <- sun_ar_auto$residuals #residuals
  res <- res - mean(res) #under the null, epsilon are mean zero
  
  ar_coef <- sun_ar_auto$model$phi #coefficients
  #Create a list with model components for arima.sim
  
  
  #Simulate series, resampling errors (should be centered)
  #and condition on p preceding values (since they are known, but not used)
  sim_unc <- function(tseries, mod, nobs, res){ 
    init <- rev(window(tseries, c(tail(time(sun),1)-nobs-mod$arma[1]+1), 
                       c(tail(time(sun),1)) - nobs)) - mod$coef["intercept"]
    mod$coef["intercept"] + filter(sample(res, size = 59, replace=TRUE), 
                                   filter = mod$model$phi, method = "recursive", sides = 1, 
                                   init = init)
  }
  
  plot(sim_unc(sun, mod = sun_ar_auto, nobs = 59, res = res), ylab = "Unconditional simulation")
  
  
  #Boostrap statistics, returned as a list here
  boot_stat <- function(tseries, npredict, p, fixorder = FALSE){
    n <- length(tseries)
    #Fit the AR model
    if(fixorder){
      ar_boot <- try(forecast::Arima(tseries[-((n-npredict+1):n)], order = c(p, 0, 0),
                                     include.mean = TRUE, include.drift = FALSE, method = "ML"))
      if(is.character(ar_boot)){
        return(list(forecast_error = rep(NA, npredict),  #forecast error
                    ar_order = p, #order of AR component
                    mu = NA #intercept
        )
        )
      }
    } else {
      ar_boot <- forecast::auto.arima(tseries[-((n-npredict+1):n)], 
                                      max.q = 0, max.Q = 0, max.d = 0, allowdrift = FALSE,
                                      max.D = 0, max.P = 0, max.p = 25,
                                      max.order = 25, stepwise = FALSE)
    }
    #Obtain forecast error for 9 periods ahead (equivalent of 1980-1988) for simulated data
    for_err <- as.vector(tseries[(n-npredict+1):n] - forecast(ar_boot, h = npredict)$mean)
    
    #Collect test statistics
    return(list(forecast_error = for_err,  #forecast error
                ar_order = ar_boot$arma[1], #order of AR component
                mu = ar_boot$coef["intercept"] #intercept
    )
    )
  }
  
  
  boot_full <- replicate(n = 199, expr = boot_stat(sim_unc(tseries = sun, 
                                                           mod = sun_ar_auto, res = res, nobs = 59), 
                                                   npredict = 9, p = p))
  
  boot_fixed <- replicate(n = 199, expr = boot_stat(sim_unc(tseries = sun, 
                                                            mod = sun_ar_auto, res = res, nobs = 59), 
                                                    npredict = 9, p = p, fixorder = TRUE))
  
  
  
  #bootstrap replicates - obtain the standard errors of the forecast errors
  for_err_full_boot <- apply(t(matrix(unlist(boot_full[1,]), nrow = 9)), 2, sd) #unconditional 
  for_err_fixed_boot <- apply(t(matrix(unlist(boot_fixed[1,]), nrow = 9)), 2, sd, na.rm = TRUE) #AR(11)
  #AR order 
  ar_order_boot <- unlist(boot_full[2,])
  plot(table(c(sun_ar_auto$arma[1], ar_order_boot))/(1+length(ar_order_boot)), ylab = "Proportion", xlab = "Order based on AIC", 
       main = "Distribution of autoregressive model order \nbased on the AIC criterion (sieve)")
  
  
  forec <- forecast(sun_ar_auto, h = 9)
  #Forecasts from forecast::forecast does not return the se
  #So reverse-engineer the calculation to retrieve those
  forec$se <- (-forec$lower[, 1] + forec$mean)/qnorm(0.5 * (1 + forec$level[1]/100))
  library(knitr)
  tab <- rbind(c(forec$se), for_err_full_boot, for_err_fixed_boot)
  row.names(tab) <- c("Nominal", "AR", "AR(11)")
  colnames(tab) <- as.character(1:9)
  kable(tab, caption = "h-step ahead prediction standard errors", digits = 2)



  
  library(boot)
  # Estimate the AR coefficients
  sun_ar <- ar(window(sun, 1930, 1979), aic = FALSE, order.max = p)
  # ar automatically selects order by AIC unless `aic = FALSE` in which case
  # it fits the model with order.max
  sun_ar$order
  [1] 11
  model <- list(ar = sun_ar$ar, order = c(p, 0, 0))
  # Statistic under study with the bootstrap Manual fitting and collection of
  # the results
  sun_fun <- function(tsb) {
    ar.fit <- ar(window(tsb, 1, 50), aic = FALSE, order.max = p)
    # Fitted using Yule-Walker equations, to avoid convergence issues and
    # because it is MUCH faster
    c(mean(tsb), c(predict(ar.fit, newdata = window(tsb, 1, 50), n.ahead = 9, 
                           se.fit = FALSE) - window(tsb, 51, 59)))
    # return prediction of time series, mean
  }
  
  # Simulation from fitted AR model, with arguments res: residuals from model
  # fit n.sim: length of series to simulate ran.args: list with components
  # `ar` and `order` From 'Bootstrap methods and their applications',
  # copyright CUP ran.gen must have precisely these arguments, in this order.
  sun_sim <- function(tseries, n.sim, ran.args) {
    rg1 <- function(n, res) {
      sample(res - mean(res), n, replace = TRUE)
    }
    ts.orig <- ran.args$ts
    ts.mod <- ran.args$model
    mean(ts.orig) + ts(arima.sim(model = ts.mod, n = n.sim, rand.gen = rg1, 
                                 res = as.vector(ran.args$res)))
  }
  # Model based bootstrap Specify the ARIMA model parameters
  sun_model <- list(order = c(sun_ar$order, 0, 0), ar = sun_ar$ar)
  sun_res <- c(scale(sun_ar$resid[!is.na(sun_ar$resid)], scale = FALSE))
  # Sieve bootstrap - also computes the test statistic on the original dataset
  # hence problems, because would usually pass residuals, and these are
  # shorter and have different time stamps use orig.t = FALSE to desactivate
  # this option
  sun_boot <- tsboot(ts(c(window(sun, 1930))), sun_fun, R = 999, sim = "model", 
                     n.sim = 59, ran.gen = sun_sim, ran.args = list(res = sun_res, ts = window(sun, 
                                                                                               1930), model = sun_model))
  
  # Standard deviations of prediction error
  apply(sun_boot$t[, -1], 2, sd)
  
  
