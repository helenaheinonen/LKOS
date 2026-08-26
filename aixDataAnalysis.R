library(copula)
library(ismev)
library(extRemes)
library(tidyverse)

source("LKOSsourceCode.R")

load("MaxMensuels.RData")
aixMonthly<-as.vector(t(maxmonth.data[18,,]))
#aixDataMonth <- data.frame(rainfall=maxmonth.data[18,])
plot(aixMonthly)

#monthlyModel.extRemes <- fevd(aixMonthly, type="GEV", method="MLE")
#summary(monthlyModel.extRemes)

monthlyModel.MLE <- gev.fit(aixMonthly) # fit GEV model to monthly maximum

#location mu: 11.707 se 0.598
#scale sigma: 10.904 se 0.501
#shape xi: 0.271 se 0.050 

# use profile LL for return level estimate 
# m = blocks for return level; since we fit monthly maxima for a k year return level we need m = k * 12 

gev.prof(monthlyModel.MLE, m = 120, xlow = 70, xup = 200, conf = 0.95, nint=1000) # 10 year return level
gev.prof(monthlyModel.MLE, m = 480, xlow = 100, xup = 300, conf = 0.95, nint=1000) # 40 year return level

gev.diag(monthlyModel.MLE) 
# diagnostic plots seem pretty good, no real reason to doubt validity of our model for just the maximum 

# qqplot curves downward at the end, so the model might be modelling a slightly heavier tail than what the data show 
# i.e. overestimate high quantiles i.e. model says rare events are more common than the data show
# makes sense though because we had a storm that would cause heavy influence (outlying point)
# this outlying point is more extreme than the model would predict 
# return level plot is contained 


## Now we'll fit largest and second largest 

load("MonthlyOrderStat.RData")
aixMonthlyOrderStat<-as.array(monthOrderStat.data[,18,,])
#station19 <- as.array(monthOrderStat.data[,19,,])
#station20 <- as.array(monthOrderStat.data[,20,,])
#station21 <- as.array(monthOrderStat.data[,21,,])

# due to seasonality changes in Rainfall we take just months Sept - Dec (9 - 12)
monthlyExtremes.septDec <- cbind(
  Year = as.vector(1976:2015),
  Largest = as.vector(aixMonthlyOrderStat[1,,9:12]),
  Second_largest = as.vector(aixMonthlyOrderStat[2,,9:12])
)


# plot all the data
ggplot(monthlyExtremes.septDec, aes(x = Year)) +
  geom_point(aes(y = monthlyExtremes.septDec[,2], colour = "Largest")) +
  geom_point(aes(y = monthlyExtremes.septDec[,3], colour = "Second largest")) +
  labs(
    x = "Year",
    y = "Rainfall",
    title = "Rainfall Extremes Sept - Dec; Aix-en-Provence",
    colour = "Monthly extremes"
  )


# fit rlarg model
septDec.model <- rlarg.fit(monthlyExtremes.septDec[,-1])
# for station 18 in Aix-en-provence
#location mu: 12.375 se 0.478
#scale sigma: 10.559 se 0.406
#shape xi: 0.406 se 0.046 
# relatively stable comapred to just the monthly maximum 
# shape grew a bit compared to 0.2 and location shifted by 1 

#rlarg.diag(septDec.model) # for some reason the QQ plot is missing -> do manually 

# Manual diagnostic plots

# QQ Plot 
par(mfrow = c(2, 2))
sorted.septDec <- apply(monthlyExtremes.septDec[,-1], 2, sort)

ss<-nrow(sorted.septDec)

q_fit1.septDec <- gevq1(septDec.model$mle, (1:ss)/(ss+1))
q_fit2.septDec <- gevq2(septDec.model$mle, (1:ss)/(ss+1))


par(mfrow = c(2, 2))
plot(
  q_fit1.septDec,
  sorted.septDec[,1],
  xlab = "Theoretical",
  ylab = "Empirical",
  main = "QQ plot (k=1)"
)
abline(0,1)

plot(
  q_fit2.septDec,
  sorted.septDec[,2],
  xlab = "Theoretical",
  ylab = "Empirical",
  main = "QQ plot (k=2)"
)
abline(0,1)

# probability plot


p_fit1.septDec <- gevf(septDec.model$mle, sorted.septDec[,1]) 
p_fit2.septDec <- gevf2(septDec.model$mle, sorted.septDec[,2])

plot(
  (1:ss)/ss,
  p_fit1.septDec,
  xlab = "Empirical",
  ylab = "Fitted",
  main = "Probability plot (k = 1)"
)
abline(0,1)


plot(
  (1:ss)/ss,
  p_fit2.septDec,
  xlab = "Empirical",
  ylab = "Fitted",
  main = "Probability plot (k = 2)"
)
abline(0,1)

# not horrible fits 

par(mfrow = c(1, 1))


empCopSample <- pobs(monthlyExtremes.septDec[,-1])
theoreticalCopSample <- rcop(ss, 2)
theoreticalCop.fin <- rcopFin(ss, k=2, blockSize = 30)


plot(empCopSample, main = "Empirical copula vs. k-extremal copula", xlab = "Largest", ylab = "Second Largest")
points(theoreticalCopSample, col="red")
legend("topleft",
       legend = c("Empirical copula", "Theoretical copula"),
       col = c("black", "red"),
       pch = 1)

modelCopSample <- cbind(
  u = gevf(monthlyOrderStat.model$mle, aixMonthlyOrderStat.sample[,2]),
  v = gevf2(monthlyOrderStat.model$mle, aixMonthlyOrderStat.sample[,3])
)

plot(empCopSample, main = "Empirical copula vs. finite block size copula", xlab = "Largest", ylab = "Second Largest")
points(theoreticalCop.fin, col = "red")
legend("topleft",
       legend = c("Empirical copula", "Theoretical copula"),
       col = c("black", "red"),
       pch = 1)
  

plot(RtransLim(empCopSample), main = "Rosenblatt Transformed Data (k-extremal copula)", xlab = "R(Largest)", ylab = "R(Second Largest)")

plot(RtransFin(empCopSample, 30), main = "Rosenblatt transformed data (finite block size copula)", xlab = "R(Largest)", ylab = "R(Second Largest)")
# not uniform 

rosenblattTest(empCopSample, finite = FALSE) # p value is 0.011
rosenblattTest(empCopSample, finite = TRUE, n = 30) # p value is 0.005

#distribution of p-values
pvals <- replicate(10, {
  rosenblattTest(empCopSample, finite = FALSE)
})

