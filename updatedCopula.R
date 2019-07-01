library(ggplot2)
library(copula)

wine_red["fixed acidity"]
wine_red[["fixed acidity"]]

wine_red_df <- as.data.frame(wine_red)
X <- as.matrix(wine_red_df)
X <- X[-c(1296, 1297), ]

mean(wine_red[["fixed acidity"]])
var(wine_red[["fixed acidity"]])

# log transform, normal dist
hist(wine_red[["fixed acidity"]])
hist(log(wine_red[["fixed acidity"]]))
X[,1] <- log(wine_red[["fixed acidity"]])

# log transform, normal dist.
hist(wine_red[["volatile acidity"]])
hist(log(wine_red[["volatile acidity"]]))
X[,2] <- log(wine_red[["volatile acidity"]])

# exponential dist.
hist((wine_red[["citric acid"]])) 

# log transform, norm dist.
hist((wine_red[["residual sugar"]]))
hist(log(wine_red[["residual sugar"]]))
X[,4] <- log(wine_red[["residual sugar"]])

hist((wine_red[["chlorides"]]))
hist(log(wine_red[["chlorides"]]))
X[,5] <- log(wine_red[["chlorides"]])

hist((wine_red[["total sulfur dioxide"]]))
hist(log(wine_red[["total sulfur dioxide"]]))
X[,6] <- log(wine_red[["total sulfur dioxide"]])

hist((wine_red[["free sulfur dioxide"]]))
hist(log(wine_red[["free sulfur dioxide"]]))
X[,7] <- log(wine_red[["free sulfur dioxide"]])

# normal dist.
hist(wine_red[["density"]])
hist(wine_red[["pH"]])

# log normal dist. 
hist(wine_red[["sulphates"]])
hist(log(wine_red[["sulphates"]]))
X[,10] <- log(wine_red[["sulphates"]])

hist((wine_red[["alcohol"]]))
hist(log(wine_red[["alcohol"]]))
X[,11] <- log(wine_red[["alcohol"]])

d <- 11 
n <- 1599

paramMargins = list(list(mean = 0, sd = 1),
                    list(mean = 0, sd = 1),
                    list(rate = 1),
                    list(mean = 0, sd = 1),
                    list(mean = 0, sd = 1),
                    list(mean = 0, sd = 1),
                    list(mean = 0, sd = 1),
                    list(mean = 0, sd = 1),
                    list(mean = 0, sd = 1),
                    list(mean = 0, sd = 1),
                    list(mean = 0, sd = 1))

# independence copula
ic <- indepCopula(dim = d)
mic <- mvdc(ic, margins = c('norm', 'norm', 'exp', 'norm', 'norm', 
                            'norm', 'norm', 'norm', 'norm', 'norm','norm'),
            paramMargins = paramMargins)
mic

start <- c(mu0 = mean(X[, 1]), sig0 = sd(X[,1]), 
           mu1 = mean(X[, 2]), sig1 = sd(X[,2]),
           lam0 = 1/mean(X[, 3]),
           mu2 = mean(X[, 4]), sig2 = sd(X[,4]),
           mu3 = mean(X[, 5]), sig3 = sd(X[,5]),
           mu4 = mean(X[, 6]), sig4 = sd(X[,6]),
           mu5 = mean(X[, 7]), sig5 = sd(X[,7]),
           mu6 = mean(X[, 8]), sig6 = sd(X[,8]),
           mu7 = mean(X[, 9]), sig7 = sd(X[,9]),
           mu8 = mean(X[, 10]), sig8 = sd(X[,10]),
           mu9 = mean(X[, 11]), sig9 = sd(X[,11]),
           th0 = 2)
start
start[is.na(start)] <- 0
start

# (mle <- fitMvdc(X[,1:11], mvdc = mic, start = start))

# Find probability based on indep. copula
Xn <- c(median(X[, 1]), median(X[, 2]), median(X[, 3]),
        median(X[, 4]), median(X[, 5]), median(X[, 6]),
        median(X[, 7]), median(X[, 8]), median(X[, 9]),
        median(X[, 10]), median(X[, 11]))

Xn

# 7th variable has NA, resulting in p7 = NA
# replace with 0
X[,7][is.na(X[,7])] <- 0

p1 <- pnorm(median(X[, 1]), mean = mean(X[, 1]), sd = sd(X[,1]), log = FALSE)
p2 <- pnorm(median(X[, 2]), mean = mean(X[, 2]), sd = sd(X[,2]), log = FALSE)
p3 <- pexp(median(X[, 3]), rate = 1/mean(X[, 3]), log = FALSE)
p4 <- pnorm(median(X[, 4]), mean = mean(X[, 4]), sd = sd(X[,4]), log = FALSE)
p5 <- pnorm(median(X[, 5]), mean = mean(X[, 5]), sd = sd(X[,5]), log = FALSE)
p6 <- pnorm(median(X[, 6]), mean = mean(X[, 6]), sd = sd(X[,6]), log = FALSE)
p7 <- pnorm(median(X[, 7]), mean = mean(X[, 7]), sd = sd(X[,7]), log = FALSE)
p8 <- pnorm(median(X[, 8]), mean = mean(X[, 8]), sd = sd(X[,8]), log = FALSE)
p9 <- pnorm(median(X[, 9]), mean = mean(X[, 9]), sd = sd(X[,9]), log = FALSE)
p10 <- pnorm(median(X[, 10]), mean = mean(X[, 10]), sd = sd(X[,10]), log = FALSE)
p11 <- pnorm(median(X[, 11]), mean = mean(X[, 11]), sd = sd(X[,11]), log = FALSE)

P <- p1 * p2 * p3 * p4 * p5 * p6 * p7 * p8 * p9 * p10 * p11
P

# clayton copula
set.seed(2019)
cc <- claytonCopula(3, dim = d)
mcc <- mvdc(cc, margins = c('norm', 'norm', 'exp', 'norm', 'norm', 
                            'norm', 'norm', 'norm', 'norm', 'norm','norm'),
            paramMargins = paramMargins)
mcc

mle <- fitMvdc(X[,1:11], mvdc = mcc, start = start)
summary(mle)

# model selection - 10-fold cross-validation
k <- 10
set.seed(2019)
xvCopula(claytonCopula(dim =d), x = X[,1:11], k = k)


# t copula
tc <- tCopula(dim = d, df = nu)
mtc <- mvdc(tc, margins = c('norm', 'norm', 'exp', 'norm', 'norm', 
                            'norm', 'norm', 'norm', 'norm', 'norm','norm'),
            paramMargins = paramMargins)
mtc

mle <- fitMvdc(X[,1:11], mvdc = mtc, start = start)
summary(mle)


# frank copula
theta <- 10 # copula parameter
fc <- frankCopula(theta, dim = d) 
set.seed(2019)

mfc <- mvdc(fc, margins = c('norm', 'norm', 'exp', 'norm', 'norm', 
                            'norm', 'norm', 'norm', 'norm', 'norm','norm'),
            paramMargins = paramMargins)

mfc

mle <- fitMvdc(X[,1:11], mvdc = mfc, start = start)
summary(mle)

# gumbel copula
gc <- gumbelCopula(2, dim = d) 
set.seed(2019)
mgc <- mvdc(gc, margins = c('norm', 'norm', 'exp', 'norm', 'norm', 
                            'norm', 'norm', 'norm', 'norm', 'norm','norm'),
            paramMargins = paramMargins)

mgc

mle <- fitMvdc(X[,1:11], mvdc = mgc, start = start)
summary(mle)

# Inference function for margins estimators

U <- cbind(pnorm(X[,1], mean = mean(X[,1]),
                 sd = sqrt((n-1)/n) * sd(X[,1])),
           pnorm(X[,2], mean = mean(X[,2]),
                 sd = sqrt((n-1)/n) * sd(X[,2])),
           pexp(X[,3], rate = 1/mean(X[,3])),
           pnorm(X[,4], mean = mean(X[,4]),
                 sd = sqrt((n-1)/n) * sd(X[,4])),
           pnorm(X[,5], mean = mean(X[,5]),
                 sd = sqrt((n-1)/n) * sd(X[,5])),
           pnorm(X[,6], mean = mean(X[,6], na.rm = TRUE),
                 sd = sqrt((n-1)/n) * sd(X[,6], na.rm = TRUE)),
           pnorm(X[,7], mean = mean(X[,7]),
                 sd = sqrt((n-1)/n) * sd(X[,7])),
           pnorm(X[,8], mean = mean(X[,8]),
                 sd = sqrt((n-1)/n) * sd(X[,8])),
           pnorm(X[,9], mean = mean(X[,9]),
                 sd = sqrt((n-1)/n) * sd(X[,9])),
           pnorm(X[,10], mean = mean(X[,10]),
                 sd = sqrt((n-1)/n) * sd(X[,10])),
           pnorm(X[,11], mean = mean(X[,11]),
                 sd = sqrt((n-1)/n) * sd(X[,11])))

# inference function bivariate RV
U <- cbind(pnorm(X[,1], mean = mean(X[,1]),
                 sd = sqrt((n-1)/n) * sd(X[,1])),
           pnorm(X[,2], mean = mean(X[,2]),
                 sd = sqrt((n-1)/n) * sd(X[,2])))


ifme <- fitCopula(claytonCopula(dim =11), data = U, start = 0.01, optim.method = 'BFGS', method = 'mpl')
summary(ifme)

# Rosenblatt tranfrom bivariate!
set.seed(2019)
cc1 <- claytonCopula(2, dim = 2)
mcc <- mvdc(cc, margins = c('norm', 'norm'),
            paramMargins = list(list(mean = 0, sd = 1),
                                list(mean = 0, sd = 1)))

start <- c(mu0 = mean(X[, 1]), sig0 = sd(X[,1]), 
           mu1 = mean(X[, 2]), sig1 = sd(X[,2]), 
           th0 = 2)

(mle <- fitMvdc(X[,1:2], mvdc = mcc, start = start))
summary(mle)

# 2.7 Rosenblatt transform
U <- rCopula(n, copula = cc1)
# apply the transformation R_C  
U. <- cCopula(U, copula = cc1)
#plot(U, xlab = quote(U*''[1]), ylab = quote(U*''[2]))
plot(U., xlab = quote(U*"'"[1]), ylab = quote(U*"'"[2]))
