################################################################################
#                         TIME SERIES ANALYSIS - HOYA 10                       #
#           Brendan Griffin, Chet Hammer, Olivia Ng, Amanda Parker             #
################################################################################

# 1. SET UP ####
library(tidyverse)
library(ggplot2)
library(forecast)
library(grid)
library(Amelia)
library(tseries)
library(scales)
library(gridExtra)
library(lmtest)
library(Rcpp)
library(questionr)
library(ggthemes)
library(Amelia)

# Load Data
liquor <- read.csv("Data/usliquorsales.csv",
                   header =TRUE)
freq.na(liquor)
sapply(liquor, class)

#Clean data
liquor$Period <- as.Date(paste(liquor$Period, "1"), 
                         "%b-%Y %d")
# Convert value to numeric
liquor$Value <- as.numeric(gsub(",","",liquor$Value))

#Converting data into a time series object
liquor_ts <-ts(liquor[,c('Value')])

# 2. PLOT DATA ####

ts_plot <- ggplot(liquor, 
                  aes(x = Period,
                      y = Value)) + 
  geom_line() + 
  labs(title = "Monthly Liquor Sales: 1992-2021", 
       x = "Year", y = "$ millions") +
  theme_igray() + 
  theme(axis.text.x = element_text(angle = 45)) + 
  scale_x_date(labels = date_format(format= "%Y"),
               breaks = date_breaks("1 year")) + 
  stat_smooth(colour = "blue") +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.title = element_text(size = 11),
        strip.background = element_rect(fill = "grey80"))
  
ts_plot

#ACF and PACF
Acf(liquor_ts,
    main = "ACF Plot")
Pacf(liquor_ts,
     main = "PACF Plot")
# ACF plots correlation coefficients between a time series and its lagged.PACF plots the partial correlation between time.
# There is a gradual decline in the ACF plot and all lags are significant. May need some differencing to achieve stationary.
# There are some significant positive and  negative correlations in the PACF plot. AR is appropriate.
# There are geometric decays in both models, which shows both seasonal and non-seasonal trend in the data.
# Lag Plot

gglagplot(liquor_ts, set.lags=1:25)


# plot 12 is highly correlated. All plots show a positive correlation.
# 3. LJUNG-BOX TEST #### 

Box.test(liquor_ts, lag=24, fitdf=0, type="Lj")

# P-value is significant which supports that ACF plot is not pure white noise 
# and there is some autocorrelation


# 4. DECOMPOSING THE TIME SERIES #### 

#Converting data into a time series object by year
liquor_tsd <-ts(liquor[,c('Value')], 
                frequency=12)

# Decomposing the time series (additive)
component.ts = decompose(liquor_tsd)

# Plot the decomposition
plot(component.ts)

# an increasing trend, seasonal variation

# Identifying and replacing outliers 
liquor$csales <-tsclean(liquor_ts)
## Do we want to remove the outliers?

#Plot the cleaned data
c_ts_plot <- ggplot(liquor, 
                    aes(x = Period,
                        y = csales)) + 
  geom_line() + 
  labs(title = "Monthly Liquor Sales: 1992-2021 (Outliers Replaced)", 
       x = "Year", y = "$ millions") +
  theme_igray() + 
  theme(axis.text.x = element_text(angle = 45)) + 
  scale_x_date(labels = date_format(format= "%Y"),
               breaks = date_breaks("1 year")) + 
  stat_smooth(colour = "blue") +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.title = element_text(size = 11),
        strip.background = element_rect(fill = "grey80"))
c_ts_plot

# compare cleaned and uncleaned plots
grid.arrange(ts_plot,
             c_ts_plot,
             ncol=1, 
             top = textGrob("Uncleaned vs Cleaned Series"))
#seem to have removed some seasonality

# Smoothing the time series decomposed time series again
my_ts <- ts(na.omit(liquor$Value), 
            frequency = 12)
plot(my_ts)

# decomposed smoothed time series
component.ts2 = decompose(my_ts)
plot(component.ts2)


# 5. FORECAST MODELS ####

# Naive Forecasting Method
naive_forecast <-naive(liquor_ts, 24)
summary(naive_forecast)
autoplot(naive_forecast)

# MAPE indicates there is 11.14% in error. 
# Check for fitted values and residuals
checkresiduals(naive_forecast)
# Residual shows homoskdasticity. Many lags indicates autocorrelation. Residual skews too the left. 

# Simple exponential smoothing
liquor_ses <- ses(my_ts, 
                  alpha=0.2, 
                  h=24)
autoplot(liquor_ses)

# Remove the trend simply by differencing the data
liquor_dif <- diff(liquor_ts)
autoplot(liquor_dif)

# Reapply the SES model
liquor_ses2 <- ses(liquor_dif, 
                   alpha=0.2, 
                   h=24)

# SES Model summary
summary(liquor_ses2)
autoplot(liquor_ses2)

# Check for fitted values and residuals
checkresiduals(liquor_ses2)

# Exponential Smoothing Model
ets_model = ets(liquor_ts, 
                allow.multiplicative.trend = TRUE)

# SES is not a good fit because there is a trend and seasonality
# ESM summary
summary(ets_model)
autoplot(ets_model)

# Check for fitted values and residuals
checkresiduals(ets_model)

# Auto ARIMA model
arima_optimal = auto.arima(my_ts,
                           seasonal = TRUE)

# ARIMA model summary
summary(arima_optimal)
autoplot(arima_optimal)

# Check for fitted values and residuals
checkresiduals(arima_optimal)

# ARIMA model with seasonal differencing
sarima <- arima(my_ts, order =c(2,1,2), 
                seasonal = list(order = c(0,1,0), period = 12))

# SARIMA model summary
summary(sarima)
autoplot(sarima)

# Check for fitted values and residuals
checkresiduals(sarima)


# 6. FORECASTING USING BEST MODEL ####

# Auto ARIMA model has the lowest AIC, lowest MAPE
# Auto ARIMA Model Forecast
forecast <- forecast(arima_optimal,
                            h = 24)
autoplot(forecast)
summary(forecast)

# Forecast values
f_values <-forecast(arima_optimal, 
                    h=24)
plot(f_values, main="")


##############################################################################################
#5. Making the series stationary (identify level of differencing required) 
#we need to remove trend by using appropriate order of difference and make the series stationary. 
#We do this by looking at acf, Dickey-Fuller Test and standard deviation.
#DICKEY FULLER TEST 
#(We have to test if Rho - 1 is significantly different than zero or not. 
#If the null hypothesis gets rejected, we’ll get a stationary time series.)
#First, confirm that the series is non-stationary using augmented DF test
adf.test(my_ts)

#To convert series to stationary, we need to know the level of differencing required
#Look at ACF (autocorrelation plot for the series to identify the order of differencing required)
Acf(my_ts)
Pacf(my_ts)

#6. Forecasting with ARIMA Model
#using differencing: lets try order 1 difference
#We will fit ARIMA(0,d,0)(0,D,0)[12] models 
#and verify acf residuals to find which ‘d’ or ‘D’ order of differencing is appropriate in our case.
#Applying only one order of difference i.e ARIMA(0,1,0)(0,0,0)
dfit1 <-arima(my_ts, order=c(0,1,0))
plot(residuals(dfit1))

Acf(residuals(dfit1))
Pacf(residuals(dfit1))

#Because the seasonal pattern is strong and stable, 
#we will want to use an order of seasonal differencing in the model. 
#Before that let’s try only with one seasonal difference i.e ARIMA(0,0,0)(0,1,0)

dfit2 <- arima(my_ts, order =c(0,0,0), seasonal = list(order = c(0,1,0), period = 12))
plot(residuals(dfit2))
Acf(residuals(dfit2))
Pacf(residuals(dfit2))
#negatives still have differencing 

#lets try and apply both seasonal and non-seasonal differencing, ARIMA(0,1,0)(0,1,0)[12]
dfit3 <- arima(my_ts, order =c(0,1,0), seasonal = list(order = c(0,1,0), period = 12))
plot(residuals(dfit3))
Acf(residuals(dfit3))
Pacf(residuals(dfit3))
#need moving average term 

#Since first ACF is -ve and most of the positive correlations are now negative (series is overdifferenced)
#we should add an MA term to the model but to know what order of MA we need,
#check the standard deviation of the models (sd=RMSE) 
summary(dfit1)
summary(dfit2)
summary(dfit3)

#We have over-differencing, so we will stop here, 
#Out of the above, dfit3 model, i.e., ARIMA(0,1,0)(0,1,0)12 has the lowest standard deviation(RMSE) and AIC. 
#Therefore, it is the correct order of differencing.
#Now, we need to identify AR/MA and SAR/SMA values and fit the model

dfit4 <- arima(my_ts, order =c(0,1,1), seasonal = list(order = c(0,1,0), period = 12))
plot(residuals(dfit4))
Acf(residuals(dfit4))
Pacf(residuals(dfit4))
#autocorrelations still not all inside

#We have over-differencing, so we will stop here, 
#Out of the above, dfit3 model, i.e., ARIMA(0,1,0)(0,1,0)12 has the lowest standard deviation(RMSE) and AIC. 
#Therefore, it is the correct order of differencing.
#Now, we need to identify AR/MA and SAR/SMA values and fit the model

#Add a one-order MA component to the seasonal part and see what we get
dfit5 <- arima(my_ts, order =c(0,1,0), seasonal = list(order = c(0,1,1), period = 12))
plot(residuals(dfit5))
Acf(residuals(dfit5))
Pacf(residuals(dfit5))

#combine a MA component to non-seasonal and one to seasonal
dfit6 <- arima(my_ts, order =c(0,1,1), seasonal = list(order = c(0,1,1), period = 12))
plot(residuals(dfit6))
Acf(residuals(dfit6))
Pacf(residuals(dfit6))

#Pending statistically significant MA coefficient and low AIC the model seems a good fit
summary(dfit4)
summary(dfit5)
summary(dfit6)

#The coeftest() function in lmtest package can help us in getting the p-values of coefficients.
coeftest(dfit6)

#Check Minimum AIC and Iterate
#We use the auto.arima() function to let R build our model with least AIC
#this function will search through combination of order parameters and provide best set
#by default it looks at maximum order of size 5 
dfit7 <- auto.arima(my_ts, seasonal = TRUE)
plot(residuals(dfit7))
Acf(residuals(dfit7))
Pacf(residuals(dfit7))
summary(dfit7)
# Smallest AIC

#7. Model Validation (n-fold holdout method)
hold <- window(ts(my_ts), start =233)

#we will forecast data for the last two years (month = 233 to 256)
fit_predicted <- arima(ts(my_ts[-c(327:351)]), order =c(2,1,1), seasonal = list(order = c(0,1,2), period = 12))

#use the model to forecast values for last 24 months. 
#Specify forecast horizon h periods ahead of prediction to be made 
#and use the fitted model to generate those predictions

forecast_pred <- forecast(fit_predicted,h=24)
plot(forecast_pred, main="")
lines(ts(my_ts))

summary(forecast_pred)

#8. Forecasting
#Next step is to forecast the sales for another 24 months ahead of time. 
f_values <-forecast(dfit7, h=24)
plot(f_values, main="")
