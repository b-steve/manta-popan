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
    ## Models with everything varying by occasion.
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

## A model.
a.model <- popanGeneral.fit.func(captlist[1:2], k=11, birthfunc = immigrationElNino.func,
                                 model=list(
                                     gp1=c("Ns.1","b1.1","b2.1",c(rep("phi1.1", 6), "phi7.1","phi7.1", rep("phi1.1", 2)), rep("p1.1", 11)),
                                     gp1=c("Ns.2","b1.1","b2.1",c(rep("phi1.2", 6), "phi7.2","phi7.2", rep("phi1.2", 2)), rep("p1.2", 11))
                                 ),
                                 ## start from the generating values: this named vector can be entered in any order
                                 startvec=c(Ns.1=1000, Ns.2=1200, b1.1=0.18, b2.1=0.1, phi1.1=0.8, phi7.1=0.4,
                                            phi1.2=0.9, phi7.2=0.5, p1.1=0.3, p1.2=0.2), printit = TRUE)

##$Ns
##     Ns.1      Ns.2 
##1234.4987  975.0703 

## A popan-General model with parameters for absolutely everything.
fit.full <- popanGeneral.fit.func(captlist[1:2], k=11, birthfunc = birth.full,
                             model=list(
                                 gp1=c("Ns.1",
                                       paste("b", 1:10, ".1", sep = ""),
                                       paste("phi", 1:10, ".1", sep = ""), 
                                       paste("p", 1:11, ".1", sep = "")),
                                 gp2=c("Ns.2",
                                       paste("b", 1:10, ".2", sep = ""),
                                       paste("phi", 1:10, ".2", sep = ""),
                                       paste("p", 1:11, ".2", sep = ""))
                             ),
                             ## start values: this named vector can be entered in any order
                             startvec=c(Ns.1=1000, Ns.2=1200,
                                        b1.1=0.1, b2.1=0.1, b3.1=0.1, b4.1=0.1, b5.1=0.1, b6.1=0.1, b7.1=0.1, b8.1=0.1, b9.1=0.1, b10.1=0.1,
                                        b1.2=0.1, b2.2=0.1, b3.2=0.1, b4.2=0.1, b5.2=0.1, b6.2=0.1, b7.2=0.1, b8.2=0.1, b9.2=0.1, b10.2=0.1,
                                        phi1.1=0.8, phi2.1=0.8, phi3.1=0.8, phi4.1=0.8, phi5.1=0.8, phi6.1=0.8, phi7.1=0.8, phi8.1=0.8, phi9.1=0.8, phi10.1=0.8,
                                        phi1.2=0.8, phi2.2=0.8, phi3.2=0.8, phi4.2=0.8, phi5.2=0.8, phi6.2=0.8, phi7.2=0.8, phi8.2=0.8, phi9.2=0.8, phi10.2=0.8,
                                        p1.1=0.3, p2.1=0.3, p3.1=0.3, p4.1=0.3, p5.1=0.3, p6.1=0.3, p7.1=0.3, p8.1=0.3, p9.1=0.3, p10.1=0.3, p11.1=0.3,
                                        p1.2=0.3, p2.2=0.3, p3.2=0.3, p4.2=0.3, p5.2=0.3, p6.2=0.3, p7.2=0.3, p8.2=0.3, p9.2=0.3, p10.2=0.3, p11.2=0.3),
                             printit = TRUE)

## A popan-General model with constant survival, because why not.
fit.2 <- popanGeneral.fit.func(captlist[1:2], k=11, birthfunc = birth.full,
                             model=list(
                                 gp1=c("Ns.1",
                                       paste("b", 1:10, ".1", sep = ""),
                                       rep("phi1.1", 10),
                                       paste("p", 1:11, ".1", sep = "")),
                                 gp2=c("Ns.2",
                                       paste("b", 1:10, ".2", sep = ""),
                                       rep("phi1.2", 10),
                                       paste("p", 1:11, ".2", sep = ""))
                             ),
                             ## start values: this named vector can be entered in any order
                             startvec=c(Ns.1=1000, Ns.2=1200,
                                        b1.1=0.1, b2.1=0.1, b3.1=0.1, b4.1=0.1, b5.1=0.1, b6.1=0.1, b7.1=0.1, b8.1=0.1, b9.1=0.1, b10.1=0.1,
                                        b1.2=0.1, b2.2=0.1, b3.2=0.1, b4.2=0.1, b5.2=0.1, b6.2=0.1, b7.2=0.1, b8.2=0.1, b9.2=0.1, b10.2=0.1,
                                        phi1.1=0.8,
                                        phi1.2=0.8,
                                        p1.1=0.3, p2.1=0.3, p3.1=0.3, p4.1=0.3, p5.1=0.3, p6.1=0.3, p7.1=0.3, p8.1=0.3, p9.1=0.3, p10.1=0.3, p11.1=0.3,
                                        p1.2=0.3, p2.2=0.3, p3.2=0.3, p4.2=0.3, p5.2=0.3, p6.2=0.3, p7.2=0.3, p8.2=0.3, p9.2=0.3, p10.2=0.3, p11.2=0.3),
                             printit = TRUE)

## Full model, with birthrates the same for males and females.
fit.3 <- popanGeneral.fit.func(captlist[1:2], k=11, birthfunc = birth.full,
                               model=list(
                                   gp1=c("Ns.1",
                                         paste("b", 1:10, ".1", sep = ""),
                                         paste("phi", 1:10, ".1", sep = ""), 
                                         paste("p", 1:11, ".1", sep = "")),
                                   gp2=c("Ns.2",
                                         paste("b", 1:10, ".1", sep = ""),
                                         paste("phi", 1:10, ".2", sep = ""),
                                         paste("p", 1:11, ".2", sep = ""))
                               ),
                             ## start values: this named vector can be entered in any order
                               startvec=c(Ns.1=1000, Ns.2=1200,
                                          b1.1=0.1, b2.1=0.1, b3.1=0.1, b4.1=0.1, b5.1=0.1, b6.1=0.1, b7.1=0.1, b8.1=0.1, b9.1=0.1, b10.1=0.1,
                                          phi1.1=0.8, phi2.1=0.8, phi3.1=0.8, phi4.1=0.8, phi5.1=0.8, phi6.1=0.8, phi7.1=0.8, phi8.1=0.8, phi9.1=0.8, phi10.1=0.8,
                                          phi1.2=0.8, phi2.2=0.8, phi3.2=0.8, phi4.2=0.8, phi5.2=0.8, phi6.2=0.8, phi7.2=0.8, phi8.2=0.8, phi9.2=0.8, phi10.2=0.8,
                                          p1.1=0.3, p2.1=0.3, p3.1=0.3, p4.1=0.3, p5.1=0.3, p6.1=0.3, p7.1=0.3, p8.1=0.3, p9.1=0.3, p10.1=0.3, p11.1=0.3,
                                          p1.2=0.3, p2.2=0.3, p3.2=0.3, p4.2=0.3, p5.2=0.3, p6.2=0.3, p7.2=0.3, p8.2=0.3, p9.2=0.3, p10.2=0.3, p11.2=0.3),
                               printit = TRUE)

## Full model, with birthrates as a function of mei.
fit.4 <- popanGeneral.fit.func(captlist[1:2], k=11, birthfunc = birth.elnino.mei,
                               model=list(
                                   gp1=c("Ns.1",
                                         paste("b", 1:2, ".1", sep = ""),
                                         paste("phi", 1:10, ".1", sep = ""), 
                                         paste("p", 1:11, ".1", sep = "")),
                                   gp2=c("Ns.2",
                                         paste("b", 1:2, ".1", sep = ""),
                                         paste("phi", 1:10, ".2", sep = ""),
                                         paste("p", 1:11, ".2", sep = ""))
                               ),
                             ## start values: this named vector can be entered in any order
                               startvec=c(Ns.1=1000, Ns.2=1200,
                                          b1.1=0.1, b2.1=0.1,
                                          phi1.1=0.8, phi2.1=0.8, phi3.1=0.8, phi4.1=0.8, phi5.1=0.8, phi6.1=0.8, phi7.1=0.8, phi8.1=0.8, phi9.1=0.8, phi10.1=0.8,
                                          phi1.2=0.8, phi2.2=0.8, phi3.2=0.8, phi4.2=0.8, phi5.2=0.8, phi6.2=0.8, phi7.2=0.8, phi8.2=0.8, phi9.2=0.8, phi10.2=0.8,
                                          p1.1=0.3, p2.1=0.3, p3.1=0.3, p4.1=0.3, p5.1=0.3, p6.1=0.3, p7.1=0.3, p8.1=0.3, p9.1=0.3, p10.1=0.3, p11.1=0.3,
                                          p1.2=0.3, p2.2=0.3, p3.2=0.3, p4.2=0.3, p5.2=0.3, p6.2=0.3, p7.2=0.3, p8.2=0.3, p9.2=0.3, p10.2=0.3, p11.2=0.3),
                               printit = TRUE)



plot(fit.full)

AIC(fit.full)
AIC(fit.2)
AIC(fit.3)

plot(fit.2)


## Try model with phi varying depending on mei.
plot(covs$mei[-1], fit.full$phis[, 1] + fit.full$rhos[, 1])

