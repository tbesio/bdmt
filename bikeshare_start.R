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

# 1.4 - Use the BH Algorithm to calculate the critical p-value for FDR of 5%
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

# 1.5 - Plot the p-value distribution
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

##### Q3: logistic regression
overload <- biketab$cnt > 500

##### Q4: treatment effects
coef(fitlin)["hum",]

##### BONUS
# Huge variety in what you could do here...













