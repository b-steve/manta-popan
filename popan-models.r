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
all.fits <- vector(mode = "list", length = n.mods)
for (i in 1:((n.mods)/100)){
    start <- 100*(i - 1) + 1
    end <- 100*(i - 1) + 100
    all.fits[start:end] <- par.fit.popan(3, arg.list = args[start:end])
    cat(i, "of", (n.mods)/100, "\n")
}
save.image("bleh.RData")
load("bleh.RData")
fits <- all.fits

probs <- which(sapply(fits, class) == "try-error")

sort(sapply(fits, AIC))[1:10]

## fits all possible models on the list above ####
#fits <- par.fit.popan(3, arg.list = args)
#save(fits, file = "2020.10.15_fits.RData")
#load("2020.10.15_fits.RData")

AIC(fits[[6290]])
plot(fits[[order(sapply(fits, AIC))[1]]])

fits[[order(sapply(fits, AIC))[4]]]$args

plot(fits[[order(sapply(fits, AIC))[1]]])

 AIC(fits[[order(sapply(fits, AIC))[4]]])

## list AICs of all possible models
sort(sapply(fits, AIC))[1:10]
## sort AICs of all models from the lowest to the highest
order(sapply(fits, AIC))[1:10]

best.fit <- do.call("fit.popan", args[[9507]])

plot(best.fit)

## best fit ####       
best.fit <- fit.popan(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                                 phi = ~ 1,
                                                                 p = ~ occasion),
                      group.pars = list(b = TRUE, phi = FALSE, p = FALSE),
                      df = covs)
AIC(best.fit)

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


