library(ggplot2)
library(copula)
library(readr)

attach(wine_red)
wine_red <- read_delim("winequality-red-1.csv", delim = ";")
wine_red["fixed acidity"]
wine_red[["fixed acidity"]]

wine_red_df <- as.data.frame(wine_red)
X. <- as.matrix(wine_red_df)

wine_red <- wine_red[-c(1296, 1297), -12]
X <- X.[-c(1296, 1297), -12]

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
X[,3] <- wine_red[["citric acid"]] + 0.000001

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

d <- ncol(X)
n <- nrow(X)

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

# Find probability based on indep. copula
Xn <- c(median(X[, 1]), median(X[, 2]), median(X[, 3]),
        median(X[, 4]), median(X[, 5]), median(X[, 6]),
        median(X[, 7]), median(X[, 8]), median(X[, 9]),
        median(X[, 10]), median(X[, 11]))

Xn

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

param_lower_bds <- lapply(mcc@margins, function(mdist) {
  if(mdist == "norm") {
    return(c(-Inf, 0))
  } else if(mdist == "exp") {
    return(0)
  }
})

param_lower_bds <- c(param_lower_bds, 0)
param_lower_bds <- unlist(param_lower_bds)

(mle <- fitMvdc(X[,1:11], mvdc = mcc, start = start))
summary(mle)

# t copula
nu <- 2
tc <- tCopula(dim = d, df = nu)
set.seed(2019)
mtc <- mvdc(tc, margins = c('norm', 'norm', 'exp', 'norm', 'norm', 
                            'norm', 'norm', 'norm', 'norm', 'norm','norm'),
            paramMargins = paramMargins)


# frank copula
theta <- 10 # copula parameter
fc <- frankCopula(theta, dim = d) 
set.seed(2019)

mfc <- mvdc(fc, margins = c('norm', 'norm', 'exp', 'norm', 'norm', 
                            'norm', 'norm', 'norm', 'norm', 'norm','norm'),
            paramMargins = paramMargins)


# gumbel copula
gc <- gumbelCopula(2, dim = d) 
set.seed(2019)
mgc <- mvdc(gc, margins = c('norm', 'norm', 'exp', 'norm', 'norm', 
                            'norm', 'norm', 'norm', 'norm', 'norm','norm'),
            paramMargins = paramMargins)

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
           pnorm(X[,6], mean = mean(X[,6]),
                 sd = sqrt((n-1)/n) * sd(X[,6])),
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

ifme <- fitCopula(claytonCopula(dim =11), data = U, start = 0.01, optim.method = 'BFGS', method = 'ml')
summary(ifme)

ifme_gumbel <- fitCopula(gumbelCopula(2, dim = 11), data = U, start = 2, optim.method = 'BFGS', method = 'ml')
summary(ifme_gumbel)

ifme_t <- fitCopula(tCopula(0.1, dim = 11, df = 2), data = U, start = c(0.1, 2), optim.method = 'BFGS', method = 'ml')
summary(ifme_t)

ifme_frank <- fitCopula(frankCopula(dim =11), data = U, start = 0.01, optim.method = 'BFGS', method = 'ml')
summary(ifme_frank)


# Maximum Pseudo-Likelihood Estimator
U. <- pobs(X[,1:11])
mple <- fitCopula(claytonCopula(dim = 11), data = U.,start = 0.01, optim.method = 'BFGS', method = "ml")
summary(mple)

mple_gumbel <- fitCopula(gumbelCopula(2, dim = 11), data = U.,start = 2, optim.method = 'BFGS', method = "ml")
summary(mple_gumbel)

mple_t <- fitCopula(tCopula(0.1, dim = 11, df = 2), data = U.,start = c(0.1, 2), optim.method = 'BFGS', method = "ml")
summary(mple_t)

mple_frank <- fitCopula(frankCopula(dim = 11), data = U.,start = 0.01, optim.method = 'BFGS', method = "ml")
summary(mple_frank)

summary(fitCopula(normalCopula(dim = 11, dispstr = "un"), data = U.))
summary(fitCopula(tCopula(dim = 11, dispstr = "un"), data = U.))

# nested Archimedean copula 
## Define a nested copula
### First nesting structure
tau1 <- c(0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 0.9) # Kendall's tau
th1 <- iTau(archmCopula('Clayton'), tau = tau1) # corresponding parameters

nlist1 <- list(th1[1], NULL, list(list(th1[2], 1, list(list(th1[3], 2:3))), # NAC structure
                                list(th1[4], 4, list(list(th1[5], 5:6))),
                                list(th1[6], 7, list(list(th1[7], 8:11)))))
NAC_1 <- onacopulaL('Clayton', nacList = nlist1) # NAC copula

### Second nesting structure
tau2 <- c(0.05, 0.1, 0.2, 0.4) # Kendall's tau
th2 <- iTau(archmCopula('Frank'), tau = tau2) # corresponding parameters

nlist2 <- list(th2[1], 1, list(list(th2[2], 2:3), # NAC structure
                             list(th2[3], 4, list(list(th2[4], 5:11)))))
NAC_2 <- onacopulaL('Frank', nacList = nlist2) # NAC copula

### Third nesting structure
tau3 <- c(0.1, 0.2, 0.5, 0.8, 0.9) # Kendall's tau
th3 <- iTau(archmCopula('Clayton'), tau = tau3) # corresponding parameters

nlist3 <- list(th3[1], 1:2, list(list(th3[2], 3:4), # NAC structure
                                list(th3[3], 5, list(list(th3[4], 6:7))), 
                                list(th3[5], 8:11)))
NAC_3 <- onacopulaL('Clayton', nacList = nlist3) # NAC copula


# IFME for nested Archimedean Copula
param <- c(.6, .3, .7, .2, .4, .8, .9)
optim_results <- optim(par = param,
                       fn = loglikNC,
                       U = U,
                       naclist = nlist1,
                       method ="BFGS",
                       control = list(fnscale = -1))

optim_results

param2 <- c(.6, .3, .7, .2)
optim_results_2 <- optim(par = param2,
                       fn = loglikNC,
                       U = U,
                       naclist = nlist2,
                       method ="BFGS",
                       control = list(fnscale = -1))

optim_results_2

# Cross Validation for nested Archimedean Copula
xvNCopula(X, nlist1, k = 10)
xvNCopula(X, nlist2, k = 10)
xvNCopula(X, nlist3, k = 10)

## Build matrix of colors
# nlist 1
cols <- matrix(1, nrow = 11, ncol = 11)
cols[1, 2:3] <- 2
cols[2, 3] <- 3
cols[4, 5:6] <- 4
cols[5:6, 5:6] <- 5
cols[7, 8:11] <- 6
cols[8:11, 8:11] <- 7
cols[lower.tri(cols)] <- t(cols)[lower.tri(cols)]
diag(cols) <- NA
splom2(U, pch = ".", pscales = 0, col.mat = cols) 

# nlist 2
cols <- matrix(1, nrow = 11, ncol = 11)
cols[2, 3] <- 2
cols[4, 5:11] <- 3
cols[5:11, 5:11] <- 4
cols[lower.tri(cols)] <- t(cols)[lower.tri(cols)]
diag(cols) <- NA
splom2(U, pch = ".", pscales = 0, col.mat = cols) 

# nlist 3
cols <- matrix(1, nrow = 11, ncol = 11)
cols[1, 2] <- 2
cols[3, 4] <- 3
cols[5, 6:7] <- 4
cols[6:7, 6:7] <- 5
cols[8:11, 8:11] <- 6
cols[lower.tri(cols)] <- t(cols)[lower.tri(cols)]
diag(cols) <- NA
splom2(U, pch = ".", pscales = 0, col.mat = cols) 


# model selection - 10-fold cross-validation
k <- 10
set.seed(2019)

xvCopula(claytonCopula(dim = d), x = X[,1:11], k = k)
xvCopula(tCopula(dim = d), x = X[,1:11], k = k)
xvCopula(frankCopula(dim = d), x = X[,1:11], k = k)
xvCopula(gumbelCopula(dim = d), x = X[,1:11], k = k)


# MLE bivariate!
set.seed(2019)
cc1 <- claytonCopula(2, dim = 2)
mcc1 <- mvdc(cc1, margins = c('norm', 'norm'),
            paramMargins = list(list(mean = 0, sd = 1),
                                list(mean = 0, sd = 1)))

start1 <- c(mu0 = mean(X[, 1]), sig0 = sd(X[,1]), 
           mu1 = mean(X[, 2]), sig1 = sd(X[,2]), 
           th0 = 2)

(mle <- fitMvdc(X[,1:2], mvdc = mcc1, start = start1))
summary(mle)

# 2.7 Rosenblatt transform
# Rosenblatt tranfrom for clayton copula
for (i in 1:11) {
  for (j in (i+1):11) {
    if(j > 11) break()
    Ui <- cbind(pnorm(X[,i], mean = mean(X[,i]),
                sd = sqrt((n-1)/n) * sd(X[,i])),
                pnorm(X[,j], mean = mean(X[,j]),
                sd = sqrt((n-1)/n) * sd(X[,j])))

    # apply the transformation R_C  
    U. <- cCopula(Ui, copula = claytonCopula(mle@estimate["th0"]), dim = 2)
    plot(U., xlab = quote(U*"'"[i]), ylab = quote(U*"'"[j]))
  }
}

# Rosenblatt tranfrom for t copula
for (i in 1:11) {
  for (j in (i+1):11) {
    if(j > 11) break()
    Ui <- cbind(pnorm(X[,i], mean = mean(X[,i]),
                      sd = sqrt((n-1)/n) * sd(X[,i])),
                pnorm(X[,j], mean = mean(X[,j]),
                      sd = sqrt((n-1)/n) * sd(X[,j])))
    
    # apply the transformation R_C  
    U. <- cCopula(Ui, copula = tCopula(mle@estimate["th0"]), dim = 2)
    plot(U., xlab = quote(U*"'"[i]), ylab = quote(U*"'"[j]))
  }
}

# Rosenblatt tranfrom for frank copula
for (i in 1:11) {
  for (j in (i+1):11) {
    if(j > 11) break()
    Ui <- cbind(pnorm(X[,i], mean = mean(X[,i]),
                      sd = sqrt((n-1)/n) * sd(X[,i])),
                pnorm(X[,j], mean = mean(X[,j]),
                      sd = sqrt((n-1)/n) * sd(X[,j])))
    
    # apply the transformation R_C  
    U. <- cCopula(Ui, copula = frankCopula(mle@estimate["th0"]), dim = 2)
    plot(U., xlab = quote(U*"'"[i]), ylab = quote(U*"'"[j]))
  }
}

# Rosenblatt tranfrom for gumbel copula
for (i in 1:11) {
  for (j in (i+1):11) {
    if(j > 11) break()
    Ui <- cbind(pnorm(X[,i], mean = mean(X[,i]),
                      sd = sqrt((n-1)/n) * sd(X[,i])),
                pnorm(X[,j], mean = mean(X[,j]),
                      sd = sqrt((n-1)/n) * sd(X[,j])))
    
    # apply the transformation R_C  
    U. <- cCopula(Ui, copula = gumbelCopula(mle@estimate["th0"]), dim = 2)
    plot(U., xlab = quote(U*"'"[i]), ylab = quote(U*"'"[j]))
  }
}

# Monte Carlo integration
mean(X[,1] <= 2)
Y <- X[,1] + X[,2] + X[,3]
mean(Y <= 1)

# Mixture Copulas
tau <- c(0.1, 0.2) # Kendall's tau
nesting_th <- iTau(archmCopula('Clayton'), tau = tau) # corresponding parameters

nesting_list1 <- list(nesting_th[1], 1:3)
nesting_list2 <- list(nesting_th[1], 1, list(list(nesting_th[2], 2:3)))
nesting_list3 <- list(nesting_th[1], 2, list(list(nesting_th[2], c(1,3))))
nesting_list4 <- list(nesting_th[1], 3, list(list(nesting_th[2], 1:2)))

param1 <- c(.6)
param2 <- c(.3, .7)
param3 <- c(.2, .4)
param4 <- c(.8, .9)

test <- list(nesting_th[1], c(6,7))
nesting_th_test <- iTau(archmCopula('Clayton'), tau = cor(U[, 6], U[, 7], method = "kendall")) # corresponding parameters
param_test <- log(0.01)
optim_results1 <- optim(par = param_test,
                       fn = loglikNC,
                       U = U,
                       naclist = test,
                       method ="BFGS",
                       control = list(fnscale = -1))

optim_results1

optim_results2 <- optim(par = param2,
                        fn = loglikNC,
                        U = U,
                        naclist = nesting_list2,
                        method ="BFGS",
                        control = list(fnscale = -1))

optim_results2

optim_results3 <- optim(par = param3,
                        fn = loglikNC,
                        U = U,
                        naclist = nesting_list3,
                        method ="BFGS",
                        control = list(fnscale = -1))

optim_results3

optim_results4 <- optim(par = param4,
                        fn = loglikNC,
                        U = U,
                        naclist = nesting_list4,
                        method ="BFGS",
                        control = list(fnscale = -1))

optim_results4

new_nesting_list1 <- updatePars( optim_results1, nesting_list1, current_index = 1)[[1]]
new_nesting_list2 <- updatePars( optim_results2, nesting_list2, current_index = 1)[[1]]
new_nesting_list3 <- updatePars( optim_results3, nesting_list3, current_index = 1)[[1]]
new_nesting_list4 <- updatePars( optim_results4, nesting_list4, current_index = 1)[[1]]

full_list <- list(new_nesting_list1, new_nesting_list2, new_nesting_list3, new_nesting_list4)
weights <- c(.1, .2, .3, .4)

optim_results <- optim(par = weights,
                       fn = loglikNCMix,
                       U = U,
                       mixlist = full_list,
                       method ="BFGS",
                       control = list(fnscale = -1))

optim_results

fitNC(10, U = U, naclist = nesting_list1)




