library(mgcv)
library(parallel)
## Sourcing model fitting functions.
source("Rfunc.R")
source("full-covs.R")
## Reading in data.
load("manta.RData")


covs$unlag.mei <- c(covs$mei[-1], 0)

covs$en.factor <- c(rep("No", 5), rep("Yes", 2), rep("No", 4))

## ## A model with fit.popan.
## fit.full <- fit.popan(captlist[1:2], model.list = list(b = ~ occasion,
##                                                        phi = ~ occasion,
##                                                        p = ~ occasion),
##                       group.pars = list(b = FALSE, phi = FALSE, p = FALSE),
##                       df = cbind(df, covs), printit = FALSE)

## Trying out par.fit.popan.
args <- list(
    ## Models with everything varying by occasion.
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ occasion,
                                                    p = ~ occasion),
         group.pars = list(b = FALSE, phi = FALSE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ occasion,
                                                    p = ~ occasion),
         group.pars = list(b = TRUE, phi = TRUE, p = TRUE),
         df = covs),
    ## Models with everything varying by occasion, but with constant
    ## survival.
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ 1,
                                                    p = ~ occasion),
         group.pars = list(b = FALSE, phi = FALSE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ 1,
                                                    p = ~ occasion),
         group.pars = list(b = TRUE, phi = TRUE, p = TRUE),
         df = covs),
    ## Modelling capture probability by survey effort.
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ 1,
                                                    p = ~ survey.eff),
         group.pars = list(b = FALSE, phi = FALSE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ 1,
                                                    p = ~ survey.eff),
         group.pars = list(b = TRUE, phi = TRUE, p = TRUE),
         df = covs),
    ## Modelling recruitment by El Nino.
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ 1,
                                                    p = ~ occasion),
         group.pars = list(b = FALSE, phi = FALSE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ 1,
                                                    p = ~ occasion),
         group.pars = list(b = TRUE, phi = TRUE, p = TRUE),
         df = covs),
    ## Modelling recruitment by lagged El Nino.
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ 1,
                                                    p = ~ occasion),
         group.pars = list(b = FALSE, phi = FALSE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ 1,
                                                    p = ~ occasion),
         group.pars = list(b = TRUE, phi = TRUE, p = TRUE),
         df = covs),
    ## Modelling capture probability by survey effort and recruitment by lagged El Nino.
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ 1,
                                                    p = ~ s(survey.eff, k = 3)),
         group.pars = list(b = FALSE, phi = FALSE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ 1,
                                                    p = ~ s(survey.eff, k = 3)),
         group.pars = list(b = TRUE, phi = TRUE, p = TRUE),
         df = covs),
    ## Final model but playing with grouping of parameters.
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ 1,
                                                    p = ~ occasion),
         group.pars = list(b = FALSE, phi = TRUE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ 1,
                                                    p = ~ occasion),
         group.pars = list(b = TRUE, phi = FALSE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ 1,
                                                    p = ~ occasion),
         group.pars = list(b = TRUE, phi = TRUE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ 1,
                                                    p = ~ occasion),
         group.pars = list(b = FALSE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ 1,
                                                    p = ~ occasion),
         group.pars = list(b = TRUE, phi = FALSE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ 1,
                                                    p = ~ occasion),
         group.pars = list(b = TRUE, phi = TRUE, p = TRUE),
         df = covs)
)

best.fit <- fit.popan(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                                 phi = ~ 1,
                                                                 p = ~ survey.eff + lag.mei),
                      group.pars = list(b = TRUE, phi = FALSE, p = FALSE),
                      df = covs)


fits <- par.fit.popan(3, arg.list = args)
sapply(fits, AIC)
order(sapply(fits, AIC))

do.call("fit.popan", args[[1]])

plot(fits[[14]])

plot(covs$unlag.mei)


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
