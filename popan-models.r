## Loading a few packages.
library(mgcv)
library(parallel)
library(R2ucare)
## Sourcing model fitting functions.
source("Rfunc.R")
source("full-covs.R")
## Reading in data.
load("manta.RData")

## Some goodness-of-fit testing.
dampier.female.hist <- dampier.captlist$female
dampier.female.freq <- rep(1, nrow(dampier.female.hist))

## Test 2 roughly tests for equal detectability: do individuals who
## were captured on occation i and are alive at occation (i + 1) have
## the same capture probability on occasion (i + 1) as individuals who
## were not captured on occasion i and are alive at occasion (i + 1)?
test2cl(dampier.female.hist, dampier.female.freq)
test2ct(dampier.female.hist, dampier.female.freq)
## Test 3 roughly tests for equal survival: are individuals who were
## first seen on occasion i seen at different rates in future than
## those individuals who were seen on occasion i, but not for the
## first time?
test3sr(dampier.female.hist, dampier.female.freq)
test3sm(dampier.female.hist, dampier.female.freq)
## An overall test incorporating all of the above.
overall_CJS(dampier.female.hist, dampier.female.freq)

dampier.male.hist <- dampier.captlist$male
dampier.male.freq <- rep(1, nrow(dampier.male.hist))

test2cl(dampier.male.hist, dampier.male.freq)
test2ct(dampier.male.hist, dampier.male.freq)
test3sr(dampier.male.hist, dampier.male.freq)
test3sm(dampier.male.hist, dampier.male.freq)
overall_CJS(dampier.male.hist, dampier.male.freq)

misool.female.hist <- misool.captlist$female
misool.female.freq <- rep(1, nrow(misool.female.hist))

test2cl(misool.female.hist, misool.female.freq)
test2ct(misool.female.hist, misool.female.freq)
test3sr(misool.female.hist, misool.female.freq)
test3sm(misool.female.hist, misool.female.freq)
overall_CJS(misool.female.hist, misool.female.freq)

misool.male.hist <- misool.captlist$male
misool.male.freq <- rep(1, nrow(misool.male.hist))

test2cl(misool.male.hist, misool.male.freq)
test2ct(misool.male.hist, misool.male.freq)
test3sr(misool.male.hist, misool.male.freq)
test3sm(misool.male.hist, misool.male.freq)
overall_CJS(misool.male.hist, misool.male.freq)

## Full Dampier model.
dampier.fit.full <- fit.popan(dampier.captlist, model.list = list(b = ~ occasion,
                                                                  phi = ~ occasion,
                                                                  p = ~ occasion),
                              group.pars = list(b = FALSE, phi = FALSE, p = FALSE),
                              group.effect = list(b = FALSE, phi = FALSE, p = FALSE))
plot(dampier.fit.full)

## Full Misool model.
misool.fit.full <- fit.popan(misool.captlist, model.list = list(b = ~ occasion,
                                                                  phi = ~ occasion,
                                                                  p = ~ occasion),
                              group.pars = list(b = FALSE, phi = FALSE, p = FALSE),
                              group.effect = list(b = FALSE, phi = FALSE, p = FALSE))
plot(misool.fit.full)

## Consider the following types of group effects for each parameter:
## - Separately estimated (group.pars = FALSE, group.effect = FALSE).
## - Group effect (group.pars = TRUE, group.effect = TRUE).
## - Same parameters (group.pars = TRUE, group.effect = FALSE).
## Consider the following models for each parameter:
## - b: ~ 1, ~ mei, ~ occasion
## - phi: ~ 1, ~ mei, ~ occasion
## - p: ~ 1, ~ mei, ~ occasion

## The object args is a list, where each component corresponds to a
## model to fit, and it contains the required arguments of fit.popan()
## for that model.
args <- list()
k <- 1
for (group.b in 1:3){
    ## Parameter grouping for recruitment.
    if (group.b == 1){
        gp.b <- FALSE
        ge.b <- FALSE
    } else if (group.b == 2){
        gp.b <- TRUE
        ge.b <- TRUE
    } else if (group.b == 3){
        gp.b <- TRUE
        ge.b <- FALSE
    }
    for (group.phi in 1:3){
        ## Parameter grouping for survival.
        if (group.phi == 1){
            gp.phi <- FALSE
            ge.phi <- FALSE
        } else if (group.phi == 2){
            gp.phi <- TRUE
            ge.phi <- TRUE
        } else if (group.phi == 3){
            gp.phi <- TRUE
            ge.phi <- FALSE
        }
        for (group.p in 1:3){
            ## Parameter grouping for detection.
            if (group.p == 1){
                gp.p <- FALSE
                ge.p <- FALSE
            } else if (group.p == 2){
                gp.p <- TRUE
                ge.p <- TRUE
            } else if (group.p == 3){
                gp.p <- TRUE
                ge.p <- FALSE
            }
            for (model.b in 1:3){
                ## Model statement for recruitment.
                if (model.b == 1){
                    m.b <- ~ 1
                } else if (model.b == 2){
                    m.b <- ~ mei
                } else if (model.b == 3){
                    m.b  <- ~ occasion
                }
                for (model.phi in 1:3){
                    ## Model statement for survival.
                    if (model.phi == 1){
                        m.phi <- ~ 1
                    } else if (model.phi == 2){
                        m.phi <- ~ mei
                    } else if (model.phi == 3){
                        m.phi  <- ~ occasion
                    }
                    for (model.p in 1:3){
                        ## Model statement for detection.
                        if (model.p == 1){
                            m.p <- ~ 1
                        } else if (model.p == 2){
                            m.p <- ~ mei
                        } else if (model.p == 3){
                            m.p  <- ~ occasion
                        }
                        args[[k]] <- list(misool.captlist,
                                          model.list = list(b = m.b,
                                                            phi = m.phi,
                                                            p = m.p),
                                          group.pars = list(b = gp.b,
                                                            phi = gp.phi,
                                                            p = gp.p),
                                          group.effect = list(b = ge.b,
                                                              phi = ge.phi,
                                                              p = ge.p),
                                          df = covs)
                        k <- k + 1        
                    }
                }
            }
        }
    }
}


fits <- par.fit.popan(3, arg.list = args)

## Print AICs of top-ten models.
sort(sapply(fits, AIC))[1:10]
## Save all AICs.
fits.AICs <- sapply(fits, AIC)

## Choose a model. m = 1 is "best" model by AIC, m = 2 is "second
## best", and so on. Try m = 1 for a model that attributes reduced
## sightings in year 8 relative to year 7 to lower survival, and m = 4
## for a model that attributes reduced sightings in year 8 relative to
## year 7 to lower detection probability.
m <- 1
fit <- fits[[order(sapply(fits, AIC))[m]]]
fit.boot <- boot.popan(fit, 1000, n.cores = 3)
plot(fit)
AIC(fit)

## Selecting fits with AIC within 10 units of the best.
best.fits <- fits[fits.AICs - min(fits.AICs) <= 10]
best.AICs <- sapply(best.fits, AIC)
## Making a dot chart for AICs.
dotchart(sort(best.AICs))
abline(v = min(best.AICs) + 10)
## Doing some bootstrap model averaging. This will take a while.
fit.ma <- boot.ma.popan(best.fits, n.cores = 3, n.boots = 1000)
## Calculating estimate summaries from the bootstrapping.
ENs.ma.weighted <- summary(fit.ma, method = "weighted")
ENs.ma.best <- summary(fit.ma, method = "best")

## P-value stuff.

## Testing total population size at year 11 vs year 1.
summary(fit.ma, method = "best", par.fun = function(x) sum(x[["ENs"]][11, ]) - sum(x[["ENs"]][1, ]),
        par.fun.p = TRUE)
## Testing average recruitment and survival over time of males vs females.
summary(fit.ma, method = "best", par.fun = function(x) diff(apply(x[["rhos"]], 2, mean)),
        par.fun.p = TRUE)





## Plotting population trajectory for the first group, with CIs, for the weighted method.
plot(ENs.ma.weighted[1:11, 1], ylim = c(0, max(ENs.ma.weighted[, 4])))
lines(ENs.ma.weighted[1:11, 1])
lines(ENs.ma.weighted[1:11, 3], lty = "dotted")
lines(ENs.ma.weighted[1:11, 4], lty = "dotted")
## ... Same for the "best" method.
points(ENs.ma.best[1:11, 1], col = "blue")
lines(ENs.ma.best[1:11, 1], col = "blue")
lines(ENs.ma.best[1:11, 3], lty = "dotted", col = "blue")
lines(ENs.ma.best[1:11, 4], lty = "dotted", col = "blue")
## ... And superimposing the bootstrap mean of the AIC weightings
## (sanity check: should be similar to the former).
points(ENs.ma.weighted[1:11, 1], col = "red")
lines(ENs.ma.weighted[1:11, 1], col = "red")

## ... And now doing the same thing for the second group.
plot(ENs.ma.weighted[12:22, 1], ylim = c(0, max(ENs.ma.weighted[, 4])))
lines(ENs.ma.weighted[12:22, 1])
lines(ENs.ma.weighted[12:22, 3], lty = "dotted")
lines(ENs.ma.weighted[12:22, 4], lty = "dotted")
points(ENs.ma.best[12:22, 1], col = "blue")
lines(ENs.ma.best[12:22, 1], col = "blue")
lines(ENs.ma.best[12:22, 3], lty = "dotted", col = "blue")
lines(ENs.ma.best[12:22, 4], lty = "dotted", col = "blue")
points(ENs.ma.weighted[12:22, 1], col = "red")
lines(ENs.ma.weighted[12:22, 1], col = "red")

## ... And now doing the same thing for total population size.
ENs.tot.ma.weighted <- summary(fit.ma, method = "weighted", pars = "ENs.tot")
ENs.tot.ma.best <- summary(fit.ma, method = "best", pars = "ENs.tot")
plot(ENs.tot.ma.weighted[, 1], ylim = c(0, max(ENs.tot.ma.weighted[, 4])))
lines(ENs.tot.ma.weighted[, 1])
lines(ENs.tot.ma.weighted[, 3], lty = "dotted")
lines(ENs.tot.ma.weighted[, 4], lty = "dotted")
points(ENs.tot.ma.best[, 1], col = "blue")
lines(ENs.tot.ma.best[, 1], col = "blue")
lines(ENs.tot.ma.best[, 3], lty = "dotted", col = "blue")
lines(ENs.tot.ma.best[, 4], lty = "dotted", col = "blue")
points(ENs.tot.ma.weighted[, 1], col = "red")
lines(ENs.tot.ma.weighted[, 1], col = "red")
