setwd("~/ArcticTemps/")
library(zoo) # data frame melting
library(qpcR) # modelling and analysis; AICC
library(dplyr) # data manipulation; might not need
library(MASS) # boxcox, confidence interval
library(forecast) # forecasting
library(GeneCycle) # periodogram
library(astsa) # times series
library(reshape2) # data frame manipulation
library(ggplot2) # prettier plots

# car: Durbin Watson Autocorrelated Errors

# data cleaning and transformation
df<-read.csv('ArcticTemps.txt', header=TRUE, sep = ",") # read in data set
sdf<-stack(df, drop=FALSE) # stack columns 
sdf$year<- seq(from = 1944, to= 2012, by=1) # create columns with years

missVal<-c(1,70,828) # mark rows with missing values (3 NAs)
sdf<-sdf[-missVal,] # remove rows with missing values

monthsAbb <- list(jan=1,feb=2,mar=3,apr=4,may=5,jun=6,jul=7,aug=8,sep=9,oct=10,nov=11,dec=12) # translate month abbrevations to numbers for easier sorting
monthNum <- function(x) {x <- tolower(x)
  sapply(x,function(x) monthsAbb[[x]])
} # function to turn month abbreviations to lower case and then replace each with it's corresponding number
sdf$ind<-monthNum(sdf$ind)

sdf<-sdf[with(sdf, order(ind)),] # order the data frame by month
sdf<-sdf[with(sdf, order(year)),] # order the data frame by year 
View(sdf)


mon.ts<- ts(sdf$values) # create the time series object
plot(y=sdf$values, x=sdf$ind, xlab="Month", ylab="Temperature", main="Average Temperatures by Month") # monthly averages
xticks <- seq(1, 12, 1)
axis(1, at = xticks, labels = xticks,  tck=1)

qqnorm(sdf$values, main="QQ-Norm Plot of Data")
qqline(sdf$values, col=2,lty=2, lwd=2)

View(sdf)

#plot(sdf[1:30,1], type='l')

# initial time series plotting
plot(mon.ts, xlab="Time", ylab="Temperature", type='l', main='Plot of Monthly Arctic Temperatures', xaxt="n")
xaxis2<- seq(1, 825, 12)
axis(1, at = xaxis2, labels = seq(1944,2012, by=1))
abline(h=mean(mon.ts), col=2, lty=2, lwd=1)
legend("bottomright", c("Monthly Data", "Overall Mean Temperature"), col=c("black", "red"), lty=1)
# Possibly strong negative autocorrelation since oscillates fairly quickly (appears choppy)
# Not stationary, there appears to be sudden spikes, as well a decrease in variability and a slight overall increase in temperature with time


####### Tests
# Augmented Dickey-Fuller Test:
adf.test(mon.ts, alternative="stationary") # p-value = 0.01
# Box-Ljung test:
Box.test(mon.ts, type="Ljung-Box") # p-value < 2.2e-16
# Box-Pierce test:
Box.test(mon.ts, type="Box-Pierce") # p-value < 2.2e-16
# p value almost zero
# stl(mon.ts) # no period
# Box.test(mon.ts, type = c("Box-Pierce", "Ljung-Box"), fitdf = 0)
shapiro.test(mon.ts)

# initial ACF and PACF
acf(mon.ts, lag.max=30, main="ACF of Original Data") # the ACF tails off very slowly; this could indicate a long memory process
pacf(mon.ts, lag.max=30, main="PACF of Original Data") # the PACF tails off somwehat after lag 18 with a few significant lags occuring after (ie: 33, 38)
# Because both the ACF and PACFs tail off slowly, the times series likely fallows an ARIMA model

# training and test sets
train.max = length(sdf$values)-40 # 40 observations is a little under 5% of the data size
m.train <- sdf[1:train.max,] # all but 40 monthly means included in training set to build model
m.test <- sdf[train.max:length(sdf$values),] # last 40 points used to test model forecasting abilities

# shifting the data up so that it is all positive 
min(m.train) # check minimum to see how much data needs to be shifted by 
m.shft <-m.train$values+21 # shift by 21 degrees so that all data is positive
int <- as.numeric(1:length(m.shft)) # time interval

######
m <-ts(m.sft)
span = 1:length(m)
bc = boxcox(m~span)
lmbda = bc$x[which(bc$y == max(bc$y))] # 2
m.bc = (m^lmbda -1 ) / lmbda
var(m.bc)
var(m)
par(mfrow = c(1,2))
plot(m.bc, main="Box Cox Transformed")
plot(m, main="Original Time Series")
par(mfrow = c(1,1))
acf(m.bc, lag.max=30, main="ACF of Box-Cox Transformed Data") 
pacf(m.bc, lag.max=30, main="PACF of Box-Cox Transformed Data")

########

m<- ts(m.shft, start=c(3,1944), frequency = 12) # turn into time series model
# box cox transformation
#fit.mon = lm(m.shft~int) # fitted model, fitting positive data by time
bc <- boxcox(m.shft~int, plotit=TRUE) # boxcox transformation (requires positive data) & plot it
# apply to find the optimal lambda so we can transform the data
lmbda = bc$x[which(bc$y == max(bc$y))]
m.bc = (bc^lmbda-1) * (1/lmbda) 
plot(m.bc)
var(m.bc)
bc^2
class(m)
int = 1:length(m)
ft = lm(m~int) # fitted model, fitting positive data by time
yval <- boxcox(m~int, plotit=TRUE) 
lmbda = yval$x[which(yval$y == max(yval$y))]
m.bc = (yval^lmbda-1) * (1/lmbda)

var(m.bc)
##########################
myModel<-arima(m.bc, order=c(8,0,1), seasonal= list(order=c(0,1,1), period=12))
myModel # gives estimated coefficients for each term
x= myModel$residuals

plot(x, main="Residuals of SARIMA (8,0,1)x(0,1,1)12 model")
shapiro.test(x) # Normality of residuals
Box.test(x, type="Ljung-Box") # Ljung Box Test: reject white noise hypothesis if p-value <0.05
Box.test(x, type="Box-Pierce") #  
Box.test((x)^2, type="Box-Pierce") # McLeod-Li: check if squares of residuals are correlated

acf(x, lag.max=20, main="ACF of Residuals")
pacf(x, lag.max=20, main="PACF of Residuals")
########################
plot(sqrt(m.shft), type='l')
plot(log(m.shft), type='l')


var(m.ts)
var(m.ts^2)
var(m.ts^.5)
var(m)

require(MASS)
m<-ts(m.train)
bcTransform <- boxcox(m~ (1:length(m))) #plots the graph
lambda = bcTransform$x[which(bcTransform$y == max(bcTransform$y))] # gives the value of ??
m.bc = (m^2-1) * (1/2)
var(m.bc)
var(m)
length(m)
length(as.numeric(m))
boxcox(mon.ts~(1:length(mon.ts), ))
m.bc = (mon.ts^2 -1 ) / 2
var(m.bc)
plot(m.ts.bc)
plot(m.ts)


lam = bc$x[which(bc$y == max(bc$y))] # x coordinate of where the curve maxes out (gives us optimal lambda)
mon.bc = (1/lam) * (m.shft^lam - 1) # transform data using lambda
op <- par(mfrow = c(1,2)) # to show both plots
ts.plot(mon.ts,main = "Original Time Series",ylab = expression(X[t]))
ts.plot(mon.bc,main = "Box-Cox Tranformed Time Series", ylab = expression(Y[t]))
par(mfrow=c(1,1)) # reset

# compare variances:
var(mon.ts) # 18.23236
var(mon.bc) # 9.332818
# decreased by 51%

# ACF and PACF of boxcox transformed data
op = par(mfrow = c(1,2))
acf(mon.ts,lag.max = 20,main = "")
pacf(mon.bc,lag.max = 50,main = "")
title("ACF/PACF of Box-Cox Transformed Arctic Data", line = -1, outer=TRUE)

# ACF and PACF of boxcox transformed data
par(mfrow=c(1,1))

# Diference at lag  = 12 (cycle determined by the ACF) to remove seasonal component
m.diff12 = diff(m, 12)
ts.plot(m.diff12, main = "Differenced at Lag 12",ylab = expression(nabla^{12}))
var(m.diff12)


# Differencing to Remove Sesonality/ Trend

m.log<- log(m.ts+22)

qqnorm(m.log, main="QQ-Norm Plot of Log Data")
qqline(m.log, col=2,lty=2, lwd=2)

m.diff1 <- diff(m.log, difference=1)
var(m.diff1)
m.ann.d <-diff(m.diff1, lag=12)
var(m.ann.d)

m.ann <-diff(m.log, lag=12)
m.mon <-diff(m.log, lag=1)
m.tri <-diff(m.log, lag=3)
m.quart <-diff(m.log, lag=4)

var(m.ts)
var(m.ann)
var(m.log)
var(m.mon)
var(m.tri)
var(m.quart)

par(mfrow = c(1,2))
plot(m.ann, main="Annually: d=12")
plot(m.mon, main="Monthlty: d=1")
plot(m.tri, main="Trimester: d=3")
plot(m.quart, main="Quarterly: d=4")



acf(m.ann, main="ACF of De-Trended TS")
pacf(m.ann, main="PACF of De-Trended TS")


ar(m.log)



m.diff <- diff(m.ts, difference=1) # difference once
plot(m.diff , ylab=expression(nabla~Y[t]), main="Differenced and Boxcox-Transformed")
abline(h=mean(m.diff), col=2, lty=2, lwd=1)
var(m.diff) # decreased from 18.66899 to 8.877436


# https://a-little-book-of-r-for-time-series.readthedocs.io/en/latest/src/timeseries.html
m.decom<- decompose(m)

m.diff2 <- diff(m.ts, difference=2) # difference once
plot(m.diff2 , ylab=expression(nabla~Y[t]), main="Differenced Data")
abline(h=mean(m.diff2), col=2, lty=2, lwd=1)
var(m.diff2) # increased from 8.877436 to 15.8395, so should not difference twice

# ACF and PACF of Box Cox Transformed, Differenced data
acf(m.diff, lag.max=100, main='ACF of Differenced Data') # still tapers off very slowly
pacf(m.diff, lag.max=150, main='PACF of Differenced Data') # still have significant lags up to 50 and at about 80, and 140

# Check variance
var(m.ts) # 18.23236
var(mon.bc) # 4324.837
var(m.diff) # 8.877436
var(m.diff2) # 15.8395
# we will use the singularly differenced data since it has the lowest variance


m.diff.12 = diff(m.ts, 12)
m.diff.12.1 <- diff(m.diff.12, difference=1) # difference once
var(m.diff.12.1)
var(m.diff.12)


# log transformation
m.log <- ts(log(m.ts+21)) # shifted by 21 to prevent negatives
plot(m.log)
# check variance 
var(m.ts) # 18.66899
var(m.log) # 0.1225263


m.log.diff <- diff(log(m.ts+21), difference=1) # take log and then difference once
var(m.log.diff) # 0.08970687

m.diff.12 = diff(log(m.ts+21), 12)
var(m.diff.12)

par(mfrow = c(1,2))
acf(m.diff12, main="ACF of Log Transformed, Diff at 12")
pacf(m.diff12, main="ACF of Log Transformed, Diff at 12")


plot(m.log.diff, main="Log and Differenced Data")
acf(m.log.diff, lag.max=100, main='ACF of Differenced Log Data') # 
pacf(m.log.diff, lag.max=150, main='PACF of Differenced Log Data') # 

m.sqrt <- sqrt(m.ts+21)
var(m.sqrt) # 0.3474679
# better than original but worse than log

# trying out various models

#### tests & diagnostic checking
# using yule walker & AIC
yw.fit <- ar(m.ts, method="yule-walker")
yw.fit #17
# using MLE & AIC
mle.fit <- ar(m.log, aic = TRUE, order.max = NULL, method="mle")
mle.fit # 12
# Shapiro 
shapiro.test(m.log)


# so either an AR value of 12 or 17
# fit = arima(y.tr, order = c(12, 1, 0), method = "ML???, xreg=1 : length(y.tr))

model1 <- arima(m.diff.12, order=c(1,0,0))
model1$aic # worst model
# attributes(model1) 
# summary(model1)
# BIC(model1)

model2 <- arima(m.diff.12, order=c(8,0,2))
model2$aic 
BIC(model2)
AICc(model2)

model3 <- arima(m.diff.12, order=c(8,0,1))
model3$aic
BIC(model3)
AICc(model3)

model4<-arima(m.diff.12, order=c(12,0,1))
model4$aic
BIC(model4)
AICc(model4)

model5<-arima(m.diff.12, order=c(17,0,1))
model5$aic
BIC(model5)
AICc(model5)

model6 <- arima(m.diff.12, order=c(8,0,2))
model6$aic
BIC(model6)
AICc(model6)

model7 <- arima(m.diff.12, order=c(25,0,1))
model7$aic
BIC(model7)
AICc(model7)

#model8 <- arima(m.diff.12, order=c(55,0,1))
#model8$aic

model9 <- arima(m.diff.12, order=c(2,0,0), seasonal = list(order = c(2, 0, 0), period = 12))
model9$aic
BIC(model9)
AICc(model9)

auto.arima(m.diff.12)

#### SARIMA
# sarima(m.ts, 8,1,1,0,1,1,12)
seasonal.ts<-arima(log(m.ts+22), order=c(8,1,2), seasonal= list(order=c(0,0,0), period=12))
seasonal.ts2<-arima(log(m.ts+22), order=c(8,1,1), seasonal= list(order=c(0,1,1), period=12))
seasonal.ts3<-arima(log(m.ts+22), order=c(5,1,1), seasonal= list(order=c(2,0,0), period=12))

seasonal.ts$aic
seasonal.ts2$aic
seasonal.ts3$aic

# arma.id(m.ts)
summary(seasonal.ts)

# Auto ARIMA
aa<-auto.arima(log(m.ts+22), trace=FALSE)
summary(aa) #ARIMA(4,1,2)  

model4<-arima(m.ts, order=c(4,1,2))
model4$aic ### best model

### Prediction

par(mfrow = c(1,1))


#m.ts
mx<-length(log(m.ts+22))
fit<-Arima(m.ts, order=c(5,1,1), seasonal= list(order=c(2,0,0), period=12), xreg=1:mx)
#fit <- Arima(m.ts, order = c(4,1,2), xreg=1:(mx-40)) # fit the arima model
m.pred <- predict(fit, n.ahead = 10, newxreg=(1+10):(mx-10)) # make predictions based off of model
upper.ci= m.pred$pred + 2*m.pred$se # find upper confidence interval
lower.ci= m.pred$pred - 2*m.pred$se # find lower confidence interval
plot(m.ts, col="black", xlim=c(160,360), ylim=c(-25,15)) # plot the original training data
lines(upper.ci, col="blue", lty="dashed") # plot upper confidence interval
lines(lower.ci, col="blue", lty="dashed") # plot lower confidence interval
points(m.pred$pred, col="red") # plot prediction points
#points((length(m.ts)+1):(length(m.ts)+40), m.pred$pred, col="red")
points((length(m.ts)+40):m.test$values, col="green") # does not show up
# autoplot(mon.ts)




m.forc<-forecast(fit, h=1, xreg=(length(m.ts)-10):(length(m.ts)))
plot(m.forc)

########## with seasonal model
sm<- seasonal.ts3
mx<-length(m.ts)
fit <- Arima(m.ts, order = c(4,1,2), xreg=1:(mx-40)) # fit the arima model
m.pred <- predict(fit, n.ahead = 40, newxreg=(1+40):(mx-40)) # make predictions based off of model
upper.ci= m.pred$pred + 2*m.pred$se # find upper confidence interval
lower.ci= m.pred$pred - 2*m.pred$se # find lower confidence interval
plot(m.ts, col="black", xlim=c(-10,850), ylim=c(-25,15)) # plot the original training data
lines(upper.ci, col="blue", lty="dashed") # plot upper confidence interval
lines(lower.ci, col="blue", lty="dashed") # plot lower confidence interval
points(m.pred$pred, col="red") # plot prediction points
#points((length(m.ts)+1):(length(m.ts)+40), m.pred$pred, col="red")
#points((length(m.ts)+40):m.test$values, col="green") # does not show up
# autoplot(mon.ts)

##############

# ar(x = diff(m.log, 1), method = "yule-walker")
fit <- Arima(m.log, order = c(5,1,1), seasonal = list(order = c(2, 0, 0), period = 12), method = "ML", xreg=1:length(m.log)) # fit the arima model
#fit2 <- Arima(m.log, order = c(25,1,1), seasonal = list(order = c(2, 0, 0), period = 12), method = "ML", xreg=1:length(m.log)) # fit the arima model
pred.tr <- predict(fit, n.ahead = 40, newxreg=(length(m.log)+1) : 230+40 ) 
U.tr= pred.tr$pred + 2*pred.tr$se
L.tr= pred.tr$pred - 2*pred.tr$se
ts.plot(m.log, xlim=c(1,230+40), ylim = c(0,max(U.tr)))
lines(U.tr, col="blue", lty="dashed")
lines(L.tr, col="blue", lty="dashed")
points((length(m.log)+1):(length(m.log)+40), pred.tr$pred, col="red")

acf(m.ann, main="ACF of De-Trended TS")
pacf(m.ann, main="PACF of De-Trended TS")



mx<-length(m.ts)
length(mon.ts)
tst <- m.test$values+22
#length(m.train)
fit <- Arima(m.ts, order = c(5,1,2), xreg=1:mx) # fit the arima model
m.pred <- predict(fit, n.ahead = 40, newxreg=(mx+1):(mx+40)) # make predictions based off of model
upper.ci= m.pred$pred + 2*m.pred$se # find upper confidence interval
lower.ci= m.pred$pred - 2*m.pred$se # find lower confidence interval
# View(mon.ts)
plot(m.ts, col="black", xlim=c(-10,850), ylim=c(-25,15)) # plot the original training data
#lines(m.ts, col="black")
lines(upper.ci, col="blue", lty="dashed") # plot upper confidence interval
lines(lower.ci, col="blue", lty="dashed") # plot lower confidence interval
lines(m.pred$pred, col="green") # plot prediction points
lines(785:825, mon.ts[785:825], col="red", pch="*")
#points((length(m.ts)+1):(length(m.ts)+40), m.pred$pred, col="red")
#points((length(m.ts)+40):m.test$values, col="green") # does not show up
# autoplot(mon.ts)