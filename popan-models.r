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
## best", and so on. Try m = 1 for a model that attributes reduced
## sightings in year 8 relative to year 7 to lower survival, and m = 4
## for a model that attributes reduced sightings in year 8 relative to
## year 7 to lower detection probability.
m <- 1
fit <- fits[[order(sapply(fits, AIC))[m]]]
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

