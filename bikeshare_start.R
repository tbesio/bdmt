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
# The estimate of sigma.hat is sqrt(SSE/(n-2))
daylm.sigma.hat <- (daylm$dev/())

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













