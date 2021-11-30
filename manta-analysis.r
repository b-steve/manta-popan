## Sourcing model fitting functions.
source("main.r")
## Reading in data.
load("misool.dampier.caphist.RData")

## Need to create an overall goodness of fit function, to do something like this:
misool.gof <- popan.gof(misool.captlist)
dampier.gof <- popan.gof(dampier.captlist)
## Getting the overall c-hat for Dampier and Misool
dampier.chat <- dampier.gof$combined$chat
misool.chat <- misool.gof$combined$chat
## You can inspect different tests like this. The tests are called
## "test2cl", "test2ct", "test3sr", and "test3sm".
## TEST2CL for group 1.
misool.gof[[1]]$test2cl
## TEST3SR for group 2.
misool.gof[[2]]$test3sr
## You can get M-arrays for each group like this:
misool.gof[[1]]$marray

fit <- fit.popan(dampier.captlist, model.list = list(b = ~ occasion, phi = ~ occasion, p = ~ occasion, ptr = ~ occasion),
                 group.pars = list(b = TRUE, phi = TRUE, p = TRUE, ptr = TRUE),
                 group.effect = list(b = FALSE, phi = FALSE, p = FALSE, ptr = FALSE),
                 df = covs)

## Doing everything for the misool analysis.
misool.wrap.out <- manta.ma.wrap(misool.captlist, mei = covs$mei, chat = misool.chat,
                                 n.boot = 5, AIC.cutoff = 10, n.cores = 3)
misool.best.fits <- misool.wrap.out$best.fits
misool.ma.fit <- misool.wrap.out$fit.ma

## Need something like this for Dampier, noting we need the chat for Dampier.
dampier.wrap.out <- manta.ma.wrap(misool.captlist, mei = covs$mei, chat = dampier.chat,
                                  n.boot = 100, AIC.cutoff = 10, n.cores = 3)
misool.best.fits <- misool.wrap.out$best.fits
misool.ma.fit <- misool.wrap.out$fit.ma

## AICs for the best Misool fits.
sapply(misool.best.fits, AIC)
## AIC for one specific model fit.
AIC(misool.best.fits[[2]])

## To see the arguments used for a specific fit you can look at the
## args component. In particular, here is the model.list component for
## one of them.

## Here's how you figure out how sex was included in the model:

## (1) If group.pars is TRUE then the coefficients for the parameter
## in model.list are shared ("grouped" across the sexes, otherwise we
## separately estimate the effects of the covariates in model.list for
## each sex.
## (2) If group.effect is TRUE, then we fit a main effect of sex. Note
## that it doesn't make sense to fit both a main effect (group.effect
## = TRUE) and also separately estimate the coefficients (group.pars =
## FALSE), so there are only three options:
## - group.pars is TRUE and group.effect is FALSE. This means males
##   and females have the same parameter values across all years.
## - group.pars is TRUE and group.effect is TRUE. This means the
##   effects of the covariates in model.list are the same for both
##   sexes, but we additionally estimate a parameter allowing a
##   constant difference across years (on the link scale) between
##   sexes.
## - group.pars is FALSE and group.effect is FALSE. This means we
##   separately fit the covariates for both sexes, so that the
##   relationship one covariate has for females might be quite
##   different than the relationship it has with males.

## Here's an example:
misool.best.fits[[2]]$args$model.list
## So we use MEI to model recruitment and survival, and have separate
## estimates per occasion for detection probabilities.
misool.best.fits[[2]]$args$group.pars
misool.best.fits[[2]]$args$group.effect

## Need something cool to happen when you do this:
## plotting Misool data
popan.plot(misool.ma.fit, year.start = 2009, year.end = 2019)
## plotting Dampier data
popan.plot(dampier.ma.fit, year.start = 2009, year.end = 2019)
## So we have:
## - Separate estimation of recruitment parameters for
##   different sexes.
## - Separate estimation of survival parameters for
##   differet sexes.
## - The detection probabilities are assumed to be the same for both
##   sexes.

## Using the shorthand notation in the paper, this model is:
## - psi(MEI*sex)phi(mei*sex)p(t).

## Makes sense when we plot it:
plot(misool.best.fits[[2]])
## - Survival and recruitment estimates have very different patterns
##   for recruitment and survival, but the pattern aligns with MEI.
## - Detection probabilities are the same for both sexes.

## A few examples of how to use the summary() function.
## Option 1: Provide a pars argument, either "ENs", "phis", "rhos",
## "ps", or "pents". Results are per-occasion, separated into the
## different groups.
summary(misool.ma.fit, pars = "ENs")
summary(misool.ma.fit, pars = "phis")
summary(misool.ma.fit, pars = "rhos")
summary(misool.ma.fit, pars = "ps")
summary(misool.ma.fit, pars = "pents")
## You can also just choose a subset of the groups.
summary(misool.ma.fit, pars = "ENs", groups = 1)
summary(misool.ma.fit, pars = "ENs", groups = 2)
## ... Or reorder them in the output for whatever reason.
summary(misool.ma.fit, pars = "ENs", groups = c(2, 1))

## Option 2: Provide a pars argument, as above, and set diffs =
## TRUE. This provides a summary of all pairwise differences between
## occasions. The output is a bit difficult to follow, however: the
## first row compares occasion 1 to occasion 1, the second row
## compares occasion 1 to occasion 2, the third row, compares occasion
## 1 to occasion 3, and so on. For 11 occasions, you will have 11*11 =
## 121 comparisons.
summary(misool.ma.fit, pars = "ENs", diffs = TRUE)

## Option 3: This is the most flexible way to use the summary()
## function. Provide a par.fun argument, which itself is a function
## that takes a fitted model and returns the function of parameters to
## summarise. The function needs to include x, and ellipses (...). You
## can access the following objects within x, each of which has a row
## for each group and a column for each occasion:
## - x$ENs
## - x$phis
## - x$rhos
## - x$ps
## - x$pents

## For example, we can replicate using pars = "ENs" for the first
## group like this:
summary(misool.ma.fit, par.fun = function(x, ...) x$ENs[, 1])
## If we just want ENs for the first group for the first two occasions:
summary(misool.ma.fit, par.fun = function(x, ...) x$ENs[1:2, 1])
## If we want the difference between the first two ENs for the first group:
summary(misool.ma.fit, par.fun = function(x, ...) x$ENs[2, 1] - x$ENs[1, 1])
## You can ask for a p-value testing the null hypothesis that the
## parameter is equal to zero by setting par.fun.p to TRUE:
summary(misool.ma.fit, par.fun = function(x, ...) x$ENs[2, 1] - x$ENs[1, 1], par.fun.p = TRUE)
## Here we're testing for a difference in ENs between groups for each occasion:
summary(misool.ma.fit, par.fun = function(x, ...) x$ENs[, 1] - x$ENs[, 2], par.fun.p = TRUE)
## Here we're getting summaries for the ratio in ENs between groups
## for each occasion:
summary(misool.ma.fit, par.fun = function(x, ...) x$ENs[, 1]/x$ENs[, 2])
## Here we're getting summaries for total estimated population sizes
## by summing over groups:
summary(misool.ma.fit, par.fun = function(x, ...) x$ENs[, 1] + x$ENs[, 2])
## Here we're testing if per-capita recruitment between occasions
## (rho) outweighs the proportion of individuals leaving the
## population between occasions (1 - phi) for the first group:
summary(misool.ma.fit, par.fun = function(x, ...) x$rhos[, 1] - (1 - x$phis[, 1]), par.fun.p = TRUE)
## And for the second group:
summary(misool.ma.fit, par.fun = function(x, ...) x$rhos[, 2] - (1 - x$phis[, 2]), par.fun.p = TRUE)
## So you can get summaries for any function of parameters you like in
## a very flexible way.
