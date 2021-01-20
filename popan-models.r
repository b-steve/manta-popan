library(mgcv)
library(parallel)
## Sourcing model fitting functions.
source("Rfunc.R")
source("full-covs.R")
## Reading in data.
load("2020.10.23_manta.RData")

enso.df <- read.csv("enso.csv")

for (i in 1:7){
    covs[, paste0("lag", i, ".mei")] <- enso.df[(8 - i):(18 - i), 2]
}



## A model with fit.popan ####
fit.full <- fit.popan(captlist[1:2], model.list = list(b = ~ 1,
                                                        phi = ~ occasion,
                                                        p = ~ occasion),
                      group.pars = list(b = TRUE, phi = TRUE, p = TRUE),
                      group.effect = list(b = TRUE, phi = TRUE, p = TRUE),
                      df = covs)

plot(fit.full)

args <- list()
k <- 1
for (group.b.i in 1:2){
    for (group.phi.i in 1:2){
        for (group.p.i in 1:2){
            for (b.i in 1:10){
                if (b.i == 1){
                    b.mod <- ~ occasion
                } else if (b.i == 2){
                    b.mod <- ~ 1
                } else if (b.i == 3){
                    b.mod <- ~ mei
                } else {
                    b.mod <- as.formula(paste0("~ ", "lag", b.i - 3, ".mei"))
                }
                for (phi.i in 1:10){
                    if (phi.i == 1){
                        phi.mod <- ~ occasion
                    } else if (phi.i == 2){
                        phi.mod <- ~ 1
                    } else if (phi.i == 3){
                        phi.mod <- ~ mei
                    } else {
                        phi.mod <- as.formula(paste0("~ ", "lag", phi.i - 3, ".mei"))
                    }
                    for (p.i in 1:19){
                        if (p.i == 1){
                            p.mod <- ~ occasion
                        } else if (p.i == 2){
                            p.mod <- ~ 1
                        } else if (p.i == 3){
                            p.mod <- ~ survey.eff
                        } else {
                            if (p.i == 4 | p.i == 5){
                                p.mei.use <- "mei"
                            } else {
                                p.mei.use <- paste0("lag", floor((p.i - 4)/2), ".mei")
                            }
                            p.mod <- as.formula(paste0("~ ", p.mei.use, " + survey.eff"[(p.i %% 2) == 1]))
                        }
                        args[[k]] <- list(caplist = captlist[1:2],
                                          model.list = list(b = b.mod,
                                                            phi = phi.mod,
                                                            p = p.mod),
                                          group.pars = list(b = TRUE,
                                                            phi = TRUE,
                                                            p = TRUE),
                                          group.effect = list(b = group.b.i == 1,
                                                              phi = group.phi.i == 1,
                                                              p = group.p.i == 1),
                                          df = covs)
                        k <- k + 1
                    }
                }
            }
        }
    }
}

## Fitting in 100-model batches.
n.mods <- length(args)
fits <- vector(mode = "list", length = n.mods)
for (i in 1:((n.mods)/100)){
    start <- 100*(i - 1) + 1
    end <- 100*(i - 1) + 100
    fits[start:end] <- par.fit.popan(3, arg.list = args[start:end])
    cat(i, "of", (n.mods)/100, "\n")
}

save.image("savepoint.RData")
## Print AICs of top-ten models.
sort(sapply(fits, AIC))[1:10]
## Save all AICs.
fits.AICs <- sapply(fits, AIC)

## Choose a model. m = 1 is "best" model by AIC, m = 2 is "second
## best", and so on.
m <- 1
fit <- fits[[order(sapply(fits, AIC))[m]]]
## Do whatever you want with this model.
fit.boot <- boot.popan(fit)
summary(fit.boot)
plot(fit)
AIC(fit)

## Playing about with model averaging using AIC weights.
av.ENs <- ma.popan(fits)
## Plotting model-averaged trajectories.
plot(av.ENs[, 1], ylim = c(0, 800))
lines(av.ENs[, 1])
points(av.ENs[, 2], pch = 2)
lines(av.ENs[, 2], lty = 2)
## Comparing to the selected model.
points(fit$ENs[, 1], col = "blue")
lines(fit$ENs[, 1], col = "blue")
points(fit$ENs[, 2], pch = 2, col = "blue")
lines(fit$ENs[, 2], lty = 2, col = "blue")

## Selecting fits with AIC within 10 units of the best.
best.fits <- fits[fits.AICs - min(fits.AICs) <= 10]
best.AICs <- sapply(best.fits, AIC)
## Making a dot chart for AICs.
dotchart(sapply(best.fits, AIC))
abline(v = min(best.AICs) + 10)
## Doing some bootstrap model averaging.
fit.ma <- ma.popan(best.fits, boot = TRUE, n.cores = 3, n.boots = 100)
## Number of times each model was chosen.
table(sapply(fit.ma, function(x) x$model.no))
## Getting averaged ENs from weighting and the bootstrap.
av.ENs.w <- ma.popan(best.fits)
av.ENs.boot <- apply(sapply(fit.ma, function(x) x$ENs), 1, mean)

## Model summary for ENs for each group from the bootstrapped model averaging.
ENs.summary <- summary(fit.ma, par.fun = function(fit) fit$ENs)
## Plotting population trajectory for the first group, with CIs.
plot(ENs.summary[1:11, 1], ylim = c(0, max(ENs.summary[, 4])))
lines(ENs.summary[1:11, 1])
lines(ENs.summary[1:11, 3], lty = "dotted")
lines(ENs.summary[1:11, 4], lty = "dotted")
## Plotting population trajectory for the second group, with CIs.
points(ENs.summary[12:22, 1], ylim = c(0, max(ENs.summary[, 4])), col = "blue")
lines(ENs.summary[12:22, 1], col = "blue")
lines(ENs.summary[12:22, 3], lty = "dotted", col = "blue")
lines(ENs.summary[12:22, 4], lty = "dotted", col = "blue")

## Model summary for the total population size from the bootstrapped model averaging.
tot.ENs.summary <- summary(fit.ma, par.fun = function(fit) apply(fit$ENs, 1, sum))
## Plotting population trajectory for the first group, with CIs.
plot(tot.ENs.summary[1:11, 1], ylim = c(0, max(tot.ENs.summary[, 4])))
lines(tot.ENs.summary[1:11, 1])
lines(tot.ENs.summary[1:11, 3], lty = "dotted")
lines(tot.ENs.summary[1:11, 4], lty = "dotted")

## Plotting model-averaged trajectories from weighted AIC.
plot(av.ENs.w[, 1], ylim = c(0, 800))
lines(av.ENs.w[, 1])
points(av.ENs.w[, 2], pch = 2)
lines(av.ENs.w[, 2], lty = 2)
## Comparing to bootstrap approach.
points(av.ENs.boot[1:11], col = "blue")
lines(av.ENs.boot[1:11], col = "blue")
points(av.ENs.boot[12:22], pch = 2, col = "blue")
lines(av.ENs.boot[12:22], lty = 2, col = "blue")
## Plotting model-averaged trajectories from weighted AIC based on all models.
points(av.ENs[, 1], ylim = c(0, 800), col = "red")
lines(av.ENs[, 1], col = "red")
points(av.ENs[, 2], pch = 2, col = "red")
lines(av.ENs[, 2], lty = 2, col = "red")
## Comparing to the best model.
points(fit$ENs[, 1], col = "green")
lines(fit$ENs[, 1], col = "green")
points(fit$ENs[, 2], pch = 2, col = "green")
lines(fit$ENs[, 2], lty = 2, col = "green")

## Bootstrapping to calculate CIs ####
## Carrying out a bootstrap procedure. This will take a while, and
## computing time will be proportional to n.boots. While you're
## playing around you can use something like n.boots = 100, but when
## you do this for realsies it's best to use 1000 or even 10000.
boot.best.fit <- boot.popan(best.fit, n.boots = 10000)

# load boot best fit
load("2020.10.23_boot-best-fit-manta.RData")

## Summary for estimated parameters.
summary.popan(boot.best.fit, function(fit) fit$fit$par)

## Summary for ENs for females.
summary.popan(boot.best.fit, function(fit) fit$ENs[, 1])

## Summary for ENs for males.
summary.popan(boot.best.fit, function(fit) fit$ENs[, 2])

## Summary for ENs for both groups.
summary.popan(boot.best.fit, function(fit) fit$ENs[, 1] + fit$ENs[, 2])



## Summary comparing first and last years (1 = females, and 2= males) #####
summary.popan(boot.best.fit, function(fit) fit$ENs[11, ] - fit$ENs[1, ])

## For p-value comparing first and last year for FEMALES:
## (1) calculate t-test statistic: estimate/std.error from the summary
372.6852/34.23032
## (2) compare to a standard normal distribution to get the p-value
## don't forget the negative value before 'estimate/std.error
## and multiply by 2 for a two-tailed test:
2*pnorm(-372.6852/34.23032)

## For p-value comparing first and last year for MALES:
## (1) calculate t-test statistic: estimate/std.error:
143.5522/59.46682
## (2) compare to a standard normal distribution to get the p-value
## and multiply by 2 for a two-tailed test:
2*pnorm(-143.5522/59.46682)

## Summary comparing females and males ####
summary.popan(boot.best.fit, function(fit) fit$ENs[, 1] - fit$ENs[, 2])

## For p-value comparing females and males in 2019 
## (if females were more abundant than males)
## (1) calculate t-test statistic: estimate/standard error
240.26506/55.57314
## (2) compare to a standard normal distribution to get the p-value
## and multiply by 2 for a two-tailed test:
2*pnorm(-240.26506/55.57314)



## Comparing males to females in terms of survival #####
summary.popan(boot.best.fit, function(fit) fit$phis[1:2, 1] - fit$phis[1:2, 2])

load("2020.10.23_manta.RData")
best.fit

# 2*pnorm(-abs(Estimate/Std Error)): use Estimate and Std Error from the summary above
# females [1,] and males [2,]
2*pnorm(-abs(0.05945483/0.02487572))
## Restricted fit with male and female survival the same.
phi.restricted.fit <- fit.popan(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                                       phi = ~ 1,
                                                                       p = ~ occasion),
                            group.pars = list(b = TRUE, phi = TRUE, p = FALSE),
                            df = covs)
## Likelihood ratio test statistics is twice difference in log-likelihoods.
phi.lrts <- 2*(phi.restricted.fit$fit$objective - best.fit$fit$objective)
phi.npar.diff <- length(best.fit$fit$par) - length(phi.restricted.fit$fit$par)
## P-value is calculated by comparing LRTS to a chi-squared distribution with npar.diff degrees of freedom.
1 - pchisq(phi.lrts, phi.npar.diff) ## p < 0.01 = females' phi is different from males' phi

## Comparing males to females in terms of detection probability
#summary.popan(boot.best.fit, function(fit) fit$ps[, 1] + fit$ps[, 2])
#2*pnorm(-abs(0.06835214/0.0192381))
## Restricted fit with male and female survival the same.
#phi.restricted.fit <- fit.popan(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
#                                                                           phi = ~ 1,
#                                                                           p = ~ occasion),
#                                group.pars = list(b = TRUE, phi = TRUE, p = FALSE),
#                                df = covs)

## Likelihood ratio test statistics is twice difference in log-likelihoods.
#lrts <- 2*(phi.restricted.fit$fit$objective - best.fit$fit$objective)
#npar.diff <- length(best.fit$fit$par) - length(phi.restricted.fit$fit$par)
## P-value is calculated by comparing LRTS to a chi-squared distribution with npar.diff degrees of freedom.
#1 - pchisq(lrts, npar.diff)


## Restricted fit with p the same for all years
p.restricted.fit <- fit.popan(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                                       phi = ~ 1,
                                                                       p = ~ 1),
                            group.pars = list(b = TRUE, phi = FALSE, p = FALSE),
                            df = covs)
## Likelihood ratio test statistics is twice difference in log-likelihoods.
p.lrts <- 2*(p.restricted.fit$fit$objective - best.fit$fit$objective)
p.npar.diff <- length(best.fit$fit$par) - length(p.restricted.fit$fit$par)
## P-value is calculated by comparing LRTS to a chi-squared distribution with npar.diff degrees of freedom.
1 - pchisq(p.lrts, p.npar.diff) ## result: 0


## Restricted fit with b the same for all years
b.restricted.fit <- fit.popan(caplist = captlist[1:2], model.list = list(b = ~ 1,
                                                                         phi = ~ 1,
                                                                         p = ~ occasion),
                              group.pars = list(b = TRUE, phi = FALSE, p = FALSE),
                              df = covs)
## Likelihood ratio test statistics is twice difference in log-likelihoods.
b.lrts <- 2*(b.restricted.fit$fit$objective - best.fit$fit$objective)
b.npar.diff <- length(best.fit$fit$par) - length(b.restricted.fit$fit$par)
## P-value is calculated by comparing LRTS to a chi-squared distribution with npar.diff degrees of freedom.
1 - pchisq(b.lrts, b.npar.diff) ## result: 0.000242 = b significantly affected by lag.mei


## Interpreting b2:
# get the values of b1 and b2 from best fits
best.fit$fit$par # use only the b2, which is b2.1 in the summary

exp(-0.4510475)

100*(1 - exp(-0.4510475))
## For every 1-unit increase in lag-mei, we estimate that the
## per-capita recruitment rate is multiplied by 0.637.

## For every 1-unit increase in lag-mei, we estimate that the
## per-capita recruitment rate decreases by 36.30%.

## save boot.best.fit to file
save(boot.best.fit, file = "2020.10.23_boot-best-fit-manta.RData")


## Plotting ENs for first group, with CIs ####
ENs.summary <- summary.popan(boot.best.fit, function(fit) fit$ENs[, 1])
## Creating a line for point estimates. The y-axis goes from 0 to the
## highest upper CI limit.
plot(ENs.summary[, 1], type = "l", ylim = c(0, max(ENs.summary[, 4])))
## Adding dotted line for lower and upper CI limits.
lines(ENs.summary[, 3], lty = "dotted")
lines(ENs.summary[, 4], lty = "dotted")


# list the model and parameter grouping ####
fits[[186]]$args$model.list
# grouping between females and males
# TRUE means FEMALES = MALES, FALSE means FEMALES are different from MALES
fits[[186]]$args$group.pars
## count the number of parameters
length(fits[[186]]$fit$par)
## AIC of the model
AIC(fits[[186]])

# plot all parameters from the best fit ####
plot(fits[[186]])

# load boot best fit data
#load("2020.09.19_boot-best-fit-manta.RData")
#load("2020.09.19_fits.RData")

## Comparing population trajectories for everything in fits.
all.ENs <- sapply(fits, function(x) apply(x$ENs, 1, sum))
plot.new()
plot.window(xlim = c(1, nrow(all.ENs)), ylim = c(0, max(all.ENs)))
box()
axis(1)
axis(2)
for (i in 1:ncol(all.ENs)){
    lines(1:nrow(all.ENs), all.ENs[, i])
}


