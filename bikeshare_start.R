### if you don't yet have data.table, run install.packages("data.table")
library(data.table)
biketab <- fread("bikeshare.csv")
# tell R which are factors
biketab[, c("dteday", "mnth","season","weekday","hr","weathersit") := list(
  factor(dteday), factor(mnth), factor(season),
  factor(weekday), factor(hr), factor(weathersit))]

### Q1: outliers and FDR
# the next command calculates total cnt by day,
# also keeping track of the corresponding yr and mnth id.
daytots <- biketab[, list(total=sum(cnt)), by=c("dteday","yr","mnth")]
row.names(daytots) <- daytots$dteday
# simple regression
daylm <- glm(total ~ yr*mnth, data=daytots)

## 1.1 Show SSE and R2
# We know that the SSE is known as the 'residual deviance' in R
# daylm$dev
# ~> [1] 726316624

# We calculate R2 as 1 - SSE/SST or 1 - residual deviance / null deviance
daylm.r2 <- 1 - daylm$dev/daylm$null.dev
# daylm.r2
# ~> [1] 0.734876

## 1.2 Code to plot the example that I use (EV for February 2012)
# plot(daytots$total[daytots$yr==1&daytots$mnth==2])
# abline(h=predict(daylm)["2012-02-01"])

## 1.3
# First calculate the standardized residuals for daylm
# The estimate of sigma.hat is sqrt(SSE/(n-k-1))
# From the summary output, we know that n-k = 707 so n-k-1 = 706
daylm.sigma.hat <- (daylm$dev/706)**0.5
# ~> [1] 1014.286

# Now we can add a column to daytots to show the standardized residuals
# First we add the predicted value to the daytots data frame
daytots$predict <- predict(daylm)
# Next we add a column for the residuals
daytots$resid <- daytots$total - daytots$predict
# Last we calculate the standardized residuals
daytots$std.resid <- daytots$resid / daylm.sigma.hat

# Next we add in the outlier p-values
daytots$pvals <- 2*pnorm(-abs(daytots$std.resid))
# And we print out the observations with the 10 smallest p-values
daytots.ordered <- daytots[order(pvals),]
head(daytots.ordered,10)

## 1.4 - Use the BH Algorithm to calculate the critical p-value for FDR of 5%
#       and list the days that are outliers based on this critical p-value
# We use the FDR utility script to find the critical p-value
source("../Utility Scripts/fdr.r")
crit.pval <- fdr_cut(daytots$pvals,.05)
# ~> [1] 2.439374e-06

# Then we print the days that would be selected as outliers using this method
daytots.outliers <- daytots.ordered[pvals<=crit.pval,]

#        dteday yr mnth total  predict     resid std.resid        pvals
# 1: 2012-10-29  1   10    22 6414.226 -6392.226 -6.302190 2.934698e-10
# 2: 2012-10-30  1   10  1096 6414.226 -5318.226 -5.243317 1.577151e-07
# 3: 2012-04-22  1    4  1027 5807.467 -4780.467 -4.713133 2.439374e-06
# So FDR would select these three days as being outliers with only a 5% false discovery rate.
# A google search showed that 10/29/2012-10/30/2012 was around the time when hurricane Sandy
# made landfall in the DC area.
# The April 22 date was less dramatic although it appears that a Nor'easter passed through the
# region on that day with potential snow/ice and cold temps.

## 1.5 - Plot the p-value distribution
# png('daytots_pvals.png')
# hist(daytots$pvals,main="Histogram of Outlier P-Values",xlab="Outlier P-Value")
# dev.off()

#### Q2: lasso regression
library(gamlr)
source("../Utility Scripts/naref.R")
mmbike <- sparse.model.matrix(
	cnt ~ . + yr*mnth + hr*notbizday,
	data=naref(biketab))[,-1]
y <- log(biketab$cnt)
## note, I need lambda.min.ratio=1e-4 because otherwise we don't get a path
## out to complex enough models (i.e. cv err is still decreasing at termination)
fitlin <- cv.gamlr( mmbike, y, lmr=1e-4, verb=TRUE )

## 2.2
# Find the OOS R2
summary(fitlin)[fitlin$seg.min,] # minimum avg OOS R2
# ~>            lambda par    oos.r2
# ~> seg69 0.001069549 669 0.9418906
summary(fitlin)[fitlin$seg.1se,] # 1se rule
# ~>            lambda par    oos.r2
# ~> seg57 0.003266248 444 0.9407428

# Plot to visualize
# png('fitlin_OOS_R2.png')
# plot(fitlin,main="OOS R2 for fitlin CV Lasso")
# dev.off()

## 2.3 - Compare the AICc, AIC, and BIC selection to each other and the CV rules
# AIC
summary(fitlin$gamlr)[which.min(AIC(fitlin$gamlr)),]
# ~>             lambda par  df        r2      aicc
# ~> seg71 0.0008879584 683 683 0.9471683 -35914.04
# AICc
summary(fitlin$gamlr)[which.min(AICc(fitlin$gamlr)),]
# ~>             lambda par  df        r2      aicc
# ~> seg71 0.0008879584 683 683 0.9471683 -35914.04
# BIC
summary(fitlin$gamlr)[which.min(BIC(fitlin$gamlr)),]
# ~>            lambda par  df        r2      aicc
# ~> seg52 0.005200791 317 317 0.9416156 -34953.39

# Probably best to show as a plot
ll <- log(fitlin$gamlr$lambda)
par(mfrow=c(1,1))
plot(fitlin$gamlr, col="grey",main="Selected Model Betas vs Log Lambda")
abline(v=ll[which.min(AICc(fitlin$gamlr))], col="black", lty=2)
abline(v=ll[which.min(AIC(fitlin$gamlr))], col="orange", lty=2)
abline(v=ll[which.min(BIC(fitlin$gamlr))], col="green", lty=2)
abline(v=log(fitlin$lambda.min), col="blue", lty=2)
abline(v=log(fitlin$lambda.1se), col="purple", lty=2)
legend("topright", bty="n", lwd=1,
  col=c("black","orange","green","blue","purple"),
  legend=c("AICc","AIC","BIC","CV.min","CV.1se"))

## 2.4 - Print the top three dteday effects
# First pull the coefficients from the 1se selected model
fitlin.1se.coef <- coef(fitlin, select="1se")
# Find the dteday coefficients and pull them into a separate structure
fitlin.dteday <- fitlin.1se.coef[grepl(glob2rx('dteday*'), rownames(fitlin.1se.coef)),]
# Store the coefficients in descending order of absolute value
fitlin.dteday.ordered <- fitlin.dteday[order(abs(fitlin.dteday),decreasing=TRUE)]
head(fitlin.dteday.ordered,3)
# ~> dteday2012-12-26 dteday2011-12-25 dteday2011-10-29
# ~>       -1.0435670       -1.0296656       -0.9707057

exp(head(fitlin.dteday.ordered,3))-1
# ~> dteday2012-12-26 dteday2011-12-25 dteday2011-10-29
# ~>       -0.6478039       -0.6428736       -0.6211844

## 2.5 - Bootstrap to the get the estimate of for the AICc and BIC selected lambdas

aicc.selected <- c()
bic.selected <- c()
n <- nrow(mmbike)

for(b in 1:100){
  ## create a matrix of resampled indices
  ib <- sample(1:n, n, replace=TRUE)
  ## create the resampled data
  mmbike.b <- mmbike[ib,]
  yb <- y[ib]

  fitb <- gamlr(mmbike.b,yb,lmr=1e-5, verb=FALSE)
  aicc.selected <- c(aicc.selected,summary(fitb)[which.min(AICc(fitb)),1])
  bic.selected <- c(bic.selected,summary(fitb)[which.min(BIC(fitb)),1])
  print(b)
}

# Plot the results
# png('bs_fitlin_aic.png')
hist(aicc.selected,main="Histogram of AICc Selected Lambdas from Bootstrap")
abline(v=summary(fitlin$gamlr)[which.min(AICc(fitlin$gamlr)),1],col='red')
# dev.off()

# png('bs_fitlin_bic.png')
hist(bic.selected,main="Histogram of BIC Selected Lambdas from Bootstrap")
abline(v=summary(fitlin$gamlr)[which.min(BIC(fitlin$gamlr)),1],col='red')
# dev.off()

##### Q3: logistic regression
# 3.1 - Creates the overload vector and converts to 0 and 1 for regression
overload <- (biketab$cnt > 500)*1

overload.fit <- gamlr(mmbike,overload,lmr=1e-4,family="binomial",standardize=FALSE)
# png('overload_lasso_path.png')
# plot(overload.fit,main="Lasso Path Plot for Overload Logistic Regression")
# dev.off()

# 3.2 - Look at the hour of day effects on probability of overload
overload.fit.coef <- coef(overload.fit)
# Find the dteday coefficients and pull them into a separate structure
overload.fit.coef.hr <- overload.fit.coef[grepl(glob2rx('hr*'), rownames(overload.fit.coef)),]

for (i in 1:length(overload.fit.coef.hr)){
  print(c(colnames(overload.fit.coef.hr)[i],overload.fit.coef.hr[i]))
}

# Focus on hr2 on a business day to compare the results to question 2
# ~> hr2 0.983
# ~> exp(0.983)-1
# ~> [1] 1.672462

## 3.3 - Find critical p cut-off

## 3.4 - Plot ROC curve
source('../Utility Scripts/roc.r')
overload.pred <- predict(overload.fit,mmbike,type="response")
# png('overload_roc_curve.png')
roc(overload.pred,overload,main="ROC Curve for Overload")

points(x= 1-mean((overload.pred<.667)[overload==0]),
y=mean((overload.pred>.667)[overload==1]),
cex=1.5, pch=20, col='red')
# dev.off()

## 3.5 - Test sample
set.seed(5807)
test <- sample(1:nrow(mmbike), 3000)

# Fit by excluding the test sample
overload.train <- overload[-test]
mmbike.train <- mmbike[-test,]
overload.fit.train <- gamlr(mmbike.train,overload.train,lmr=1e-4,family="binomial",standardize=FALSE)

# Predict for the test set
overload.pred.test <- predict(overload.fit.train,mmbike[test,],type="response")

# Plot ROC curve for OOS predictions
# png('oos_roc.png')
# roc(overload.pred.test,overload[test],main="ROC Curve for OOS Prediction of Overload")
# dev.off()

##### Q4: treatment effects
## 4.1 - Impact of humidity baised on 'naive' regression
coef(fitlin)["hum",]
# ~> [1] -0.04980466
exp(coef(fitlin)["hum",])-1
# ~> [1] -0.04858474

## 4.2 - Predict humidity from the other covariates
d <- biketab[,hum]
covars <- subset(biketab, select = c(1:12))
controls <- sparse.model.matrix(hum ~ .,
  data=naref(covars))[,-1]

# Run a regression to predict d from the controls
treat <- gamlr(controls,d)

# Predict dhat from the fitted betas
dhat <- predict(treat, controls, type="response")
plot(dhat,d,bty="n",pch=21,bg=8)

# Calculate the in sample R2 for the d on x regression
cor(drop(dhat),d)^2
# ~> [1] 0.8107548

## 4.3 - Estimate independent effect of humidity
causal <- gamlr(cBind(d,dhat,controls),y,free=2)

# What is the beta on the independent effect of d?
coef(causal)["d",]
# ~> [1] -0.03802219

## 4.4 Extend the fitlin to include a temperature / humidity interaction
mmbike.th <- sparse.model.matrix(
  cnt ~ . + yr*mnth + hr*notbizday + temp*hum,
  data=naref(biketab))[,-1]

fitlin.th <- cv.gamlr( mmbike.th, y, lmr=1e-4, verb=TRUE )

coef(fitlin.th$gamlr)["hum",]
# ~> [1] -0.04803867
coef(fitlin.th$gamlr)["temp:hum",]
# ~> [1] 0.05051735

## 4.5 - Extend 4.3 to include the temp hum interaction
causal.dep <- cbind(d,dhat,covars)
causal.dep.int <- sparse.model.matrix(d ~ . + temp*dhat - dhat,data=naref(causal.dep))

causal.temp <- gamlr(causal.dep.int,y,free=)

##### BONUS
# Huge variety in what you could do here...













