library(MASS)
library(leaps)
library(glmnet)

# load data
Data <- read.table("UsedCar2.csv",sep=",",header=FALSE)
names(Data) <- c("Y","X1","X2","X3","X4","X5","X6")

# simply do linear regression on all variables
fit <- lm(Y~X1+X2+X3+X4+X5+X6,data=Data)
summary(fit)
cor(Data)

# draw residual plot
plot(fit$fitted.values,fit$residuals)
plot(Data$X1,fit$residuals)
plot(Data$X2,fit$residuals)
plot(Data$X3,fit$residuals)
plot(Data$X4,fit$residuals)
plot(Data$X5,fit$residuals)
plot(Data$X6,fit$residuals)
# see there are several outliers

# test if outliers are influential use cook's distance
outliers <- which(cooks.distance(fit)>0.1)
cooks.distance(fit)[c(138,139)]
Data[c(138,139),]
# find case 138 & 139 are influential outliers
#( Here should explain why they are outliers, recording error etc.)

# remove case 138 & 139 from data 
DataWithoutOutliers <- Data[-c(138,139),]

# use data without outliers to do the linear regression 
fit2 <- lm(Y~X1+X2+X3+X4+X5+X6,data=DataWithoutOutliers)
summary(fit2)
plot(fit2$fitted.values,fit2$residuals)
plot(DataWithoutOutliers$X1,fit2$residuals)
plot(DataWithoutOutliers$X2,fit2$residuals)
plot(DataWithoutOutliers$X3,fit2$residuals)
plot(DataWithoutOutliers$X4,fit2$residuals)
plot(DataWithoutOutliers$X5,fit2$residuals)
plot(DataWithoutOutliers$X6,fit2$residuals)
X1 <- DataWithoutOutliers$X1
X2 <- DataWithoutOutliers$X2
X3 <- DataWithoutOutliers$X3
X4 <- DataWithoutOutliers$X4
X5 <- DataWithoutOutliers$X5
X6 <- DataWithoutOutliers$X6
Y <- DataWithoutOutliers$Y
n <- dim(DataWithoutOutliers)[1]
# from the plot we find error with unconstant variance

###############################
# do Brown-Forsythe Tests to test if errors have unconstant variance
# here we do tests with all predictor variables
v = c(2:7)
for (i in v){
  #1. Break the residuals into two groups. 
  Group1 <- fit2$residuals[DataWithoutOutliers[,i]<mean(DataWithoutOutliers[,i])]
  Group2 <- fit2$residuals[DataWithoutOutliers[,i]>=mean(DataWithoutOutliers[,i])]
  #2. Obtain the median of each group, using the commands: 
  M1 <- median(Group1) 
  M2 <- median(Group2) 
  #3. Obtain the mean absolute deviation for each group, using the commands: 
  D1 <- sum( abs( Group1 - M1 )) / length(Group1) 
  D2 <- sum( abs( Group2 - M2 )) / length(Group2)
  #4. Calculate the pooled standard error, using the command: 
  s <- sqrt( ( sum( ( abs(Group1 - M1) - D1 )^2 ) + sum( ( abs(Group2 - M2) - D2 )^2 ) ) / (n-2) ) 
  #5. Finally, calculate the Brown-Forsythe test statistic, using the command: 
  t.bf <- ( D1 - D2 ) / ( s * sqrt( 1/length(Group1) + 1/length(Group2) ) ) 
  #output test statistic
  print(c(i,t.bf))
}
#6 Once you obtain this value, you can compare it to the critical value for any given alpha level to determine whether or not to conclude constancy of error variance, 
# or you can find its P-value. 
alpha <- 0.05
qt(1-alpha/2, n-2)   # find the catical value
# And the P-value can be found by typing: 
2*(1-pt( abs(t.bf), n-2))
# find error with X3 has non-constant variance

# do the box-cox transformation
boxcox(Y~X1+X2+X3+X4+X5+X6)

# find lambda close to 0, seems like 0.25

# use lamda=0 , i.e. log
fit3 <- lm(log(Y)~X1+X2+X3+X4+X5+X6,data=DataWithoutOutliers)
summary(fit3)
plot(fit3$fitted.values,fit3$residuals)
plot(DataWithoutOutliers$X1,fit3$residuals)
plot(DataWithoutOutliers$X2,fit3$residuals)
anova(fit3)

# use lamda=0.25
Y.boxcox <- (Y^0.2-1)/0.25
fit4 <- lm(Y.boxcox~X1+X2+X3+X4+X5+X6,data=DataWithoutOutliers)
summary(fit4)
#find log transformation has larger R-squared, so use fit3 instead of fit4

#################################
# see correlation matrix
cor(DataWithoutOutliers)
# find X1&X2, X1&X5,X2&X5 have high correlation

# calculate VIF of all predictor variables
VIF.fun <- function(Y, X, DataWithoutOutliers)
{
  fit.vif <- lm(Data[,Y] ~ Data[,X[1]] + Data[,X[2]]+Data[,X[3]]+Data[,X[4]]+Data[,X[5]])
  fit.vif.r2 <- summary(fit.vif)$r.squared
  VIF <- 1/(1-fit.vif.r2)
  VIF
}
VIF <- rep(0,6)
VIF[1] <- VIF.fun(2, c(1,3,4,5,6,7), DataWithoutOutliers)
VIF[2] <- VIF.fun(3, c(1,2,4,5,6,7), DataWithoutOutliers)
VIF[3] <- VIF.fun(4, c(1,2,3,5,6,7), DataWithoutOutliers)
VIF[4] <- VIF.fun(5, c(1,2,3,4,6,7), DataWithoutOutliers)
VIF[5] <- VIF.fun(6, c(1,2,3,4,5,7), DataWithoutOutliers)
VIF[6] <- VIF.fun(7, c(1,2,3,4,5,6), DataWithoutOutliers)
VIF
# find X2 and X3 has larger VIF, may have multicollinearity 
# (do Ridge regression at last)

###################################
# X3 and X6 has large p-value in fit3
# test if X3 can be dropped ( beta3 = 0)
alpha <- 0.05
n <- dim(DataWithoutOutliers)[1]
qt(1-alpha/2,n-2)

fit5 <- lm(log(Y)~X1+X2+X4+X5+X6,data=DataWithoutOutliers)
anova(fit5,fit3)
summary(fit5)
# X3 can be dropped

# test if X6 can be dropped ( beta4 = beta5 = beta6 = 0 )
fit6 <- lm(log(Y)~X1+X2,data=DataWithoutOutliers)
anova(fit6,fit5)
# X4,X5,X6 cannot be dropped

###############################
# add second-order term and see if they can be dropped
fit7 <- lm(log(Y)~X1+X2+X4+X5+X6+X1*X2+X1*X3+X1*X4+X1*X5+X1*X6
           +X2*X3+X2*X4+X2*X5+X2*X6+X3*X4+X3*X5+X3*X6+
             I(X1^2)+I(X2^2)+I(X3^2)+I(X4^2)+I(X5^2)+I(X6^2),data=DataWithoutOutliers)
anova(fit7,fit5)
summary(fit7)
# second-order term cannot be dropped

# model selection

# Cp criterion
leaps( x=DataWithoutOutliers[,2:7], y=log(DataWithoutOutliers[,1]),
       names=names(DataWithoutOutliers)[2:7],method="Cp")
# the best is without X3&X6, without X3

# adjust R-squared criterion
Result.adjr2 <- leaps( x=DataWithoutOutliers[,2:7], y=log(DataWithoutOutliers[,1]),
       names=names(DataWithoutOutliers)[2:7],method="adjr2")
max(Result.adjr2$adjr2)
# the best is Model37. Model 27 is also a good choice

# stepwise selection
Base <- lm( Y ~ 1, data=DataWithoutOutliers )
step(Base, scope = list( upper=fit7, lower=~1 ), direction = "both", trace=TRUE)

# add X2^2,X1*X4 based on fit5
fit8 <- lm(log(Y)~X1+X2+X4+X5+X6+I(X2^2)+X1*X4,data=DataWithoutOutliers)
summary(fit8)

##################################
# compare fit5 & fit8 with ridge regression respectively

#transform of Y
DataWithoutOutliers$logY = log(DataWithoutOutliers$Y)

# Set random seed
set.seed(1)

# split whole data into train data and test data 
train = runif(nrow(DataWithoutOutliers)) < .5

# freq of train
table(train)

################################
# fit your start model - fit5
fit9 = lm(data = DataWithoutOutliers, subset = train, logY ~ . -X3 - Y)
summary(fit9)

# predict on your test data
yhat_1 = predict(fit9, DataWithoutOutliers[!train,])

# get the MSE
testmse.fit9 <- mean((DataWithoutOutliers$logY[!train] - yhat_1)^2)

# get the residual plot 
plot(fit9)

#########################
# create model matrix
x = model.matrix(logY ~ . - X3 - Y, data = DataWithoutOutliers)

# fit ridge regression
fit.ridge = glmnet(x[train,], DataWithoutOutliers$logY[train], alpha = 0)

# ridge trace
plot(fit.ridge, xvar = "lambda")

# use cv to select best ridge model
fit.cv = cv.glmnet(x[train,], DataWithoutOutliers$logY[train], alpha = 0)
fit.cv$lambda.min

# get test MSE
yhat_3 = predict(fit.ridge, s = fit.cv$lambda.min, newx = x[!train,])
testmse.fit5.ridge <- mean((DataWithoutOutliers$logY[!train] - yhat_3)^2)


coef(fit.cv,s="lambda.min")
##################################
# fit your start model - fit8
fit10 = lm(data = DataWithoutOutliers, subset=train, logY ~ . -X3 - Y +I(X2^2)+X1*X4)
summary(fit10)

# predict on your test data
yhat_1 = predict(fit10, DataWithoutOutliers[!train,])

# get the MSE
testmse.fit10 <- mean((DataWithoutOutliers$logY[!train] - yhat_1)^2)

# get the residual plot 
plot(fit10)

# create model matrix
z = model.matrix(logY ~ . - X3 - Y + I(X2^2) + X1*X4, data = DataWithoutOutliers)

# fit ridge regression
fit.ridge = glmnet(z[train,], DataWithoutOutliers$logY[train], alpha = 0)

# ridge trace
plot(fit.ridge, xvar = "lambda")

# use cv to select best ridge model
fit.cv = cv.glmnet(z[train,], DataWithoutOutliers$logY[train], alpha = 0)
fit.cv$lambda.min

# get test MSE
yhat_3 = predict(fit.ridge, s = fit.cv$lambda.min, newx = z[!train,])
testmse.fit8.ridge <- mean((DataWithoutOutliers$logY[!train] - yhat_3)^2)

c(testmse.fit9,testmse.fit5.ridge)
c(testmse.fit10,testmse.fit8.ridge)


# find fit5.ridge will be the best one.

# interpretation  : equation, multicollinearity, overall significance&individual significance
# variable interpretation,commonsense question
