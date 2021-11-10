## Sourcing model fitting functions.
source("main.r")
## Reading in data.
load("manta.RData")

## Need to create an overall goodness of fit function, to do something like this:
## misool.gof <- popan.gof(misool.captlist)
## dampier.gof <- popan.gof(dampier.captilst)
## dampier.chat <- dampier.gof$chat

## Doing everything for the misool analysis.
misool.ma.fit <- manta.ma.wrap(misool.captlist, mei = covs$mei, chat = 1,
                               n.boot = 100, AIC.cutoff = 10, n.cores = 3)
## Need something like this for Dampier, noting we need the chat for Dampier.
## dampier.ma.fit <- manta.ma.wrap(misool.captlist, mei = covs$mei, chat = dampier.chat,
##                                 n.boot = 100, AIC.cutoff = 10, n.cores = 3)

## Need something cool to happen when you do this:
## plot(misool.fit)

## Need the output for the summary() function tidied up a bit, to show
## more clearly what the estimates are. At the moment rows 1-11 are
## estimated abundances for females, 2-12 are for males, but it's not
## clear.
summary(misool.ma.fit)
