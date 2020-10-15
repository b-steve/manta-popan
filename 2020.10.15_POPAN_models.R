library(mgcv)
library(parallel)
## Sourcing model fitting functions.
source("Rfunc.R")
source("full-covs.R")
## Reading in data.
load("2020.10.10_manta.RData")

covs$unlag.mei <- c(covs$mei[-1], 0)
covs$en.factor <- c(rep("No", 5), rep("Yes", 2), rep("No", 4))

## A model with fit.popan ####
#fit.full <- fit.popan(captlist[1:2], model.list = list(b = ~ occasion,
#                                                        phi = ~ occasion,
#                                                        p = ~ occasion),
#                       group.pars = list(b = FALSE, phi = FALSE, p = FALSE),
#                       df = cbind(df, covs), printit = FALSE)

## Trying out par.fit.popan ####
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
         group.pars = list(b = FALSE, phi = FALSE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ occasion,
                                                    p = ~ occasion),
         group.pars = list(b = FALSE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ occasion,
                                                    p = ~ occasion),
         group.pars = list(b = TRUE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ occasion,
                                                    p = ~ occasion),
         group.pars = list(b = TRUE, phi = TRUE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ occasion,
                                                    p = ~ occasion),
         group.pars = list(b = TRUE, phi = FALSE, p = FALSE),
         df = covs),
    
    ## Models with everything varying by occasion,
    ## except for p (by survey effort)
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ occasion,
                                                    p = ~ survey.eff),
         group.pars = list(b = FALSE, phi = FALSE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ occasion,
                                                    p = ~ survey.eff),
         group.pars = list(b = FALSE, phi = FALSE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ occasion,
                                                    p = ~ survey.eff),
         group.pars = list(b = FALSE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ occasion,
                                                    p = ~ survey.eff),
         group.pars = list(b = TRUE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ occasion,
                                                    p = ~ survey.eff),
         group.pars = list(b = TRUE, phi = TRUE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ occasion,
                                                    p = ~ survey.eff),
         group.pars = list(b = TRUE, phi = FALSE, p = FALSE),
         df = covs),
    
    ## Models with everything varying by occasion, except for p (by MEI)
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ occasion,
                                                    p = ~ mei),
         group.pars = list(b = FALSE, phi = FALSE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ occasion,
                                                    p = ~ mei),
         group.pars = list(b = FALSE, phi = FALSE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ occasion,
                                                    p = ~ mei),
         group.pars = list(b = FALSE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ occasion,
                                                    p = ~ mei),
         group.pars = list(b = TRUE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ occasion,
                                                    p = ~ mei),
         group.pars = list(b = TRUE, phi = TRUE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ occasion,
                                                    p = ~ mei),
         group.pars = list(b = TRUE, phi = FALSE, p = FALSE),
         df = covs),
    
    ## Models with everything varying by occasion, except for p (by lag.mei)
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ occasion,
                                                    p = ~ lag.mei),
         group.pars = list(b = FALSE, phi = FALSE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ occasion,
                                                    p = ~ lag.mei),
         group.pars = list(b = FALSE, phi = FALSE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ occasion,
                                                    p = ~ lag.mei),
         group.pars = list(b = FALSE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ occasion,
                                                    p = ~ lag.mei),
         group.pars = list(b = TRUE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ occasion,
                                                    p = ~ lag.mei),
         group.pars = list(b = TRUE, phi = TRUE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ occasion,
                                                    p = ~ lag.mei),
         group.pars = list(b = TRUE, phi = FALSE, p = FALSE),
         df = covs),
    
    ## Models with everything varying by occasion
    ## except for p (by survey effort + mei)
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ occasion,
                                                    p = ~ survey.eff + mei),
         group.pars = list(b = FALSE, phi = FALSE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ occasion,
                                                    p = ~ survey.eff + mei),
         group.pars = list(b = FALSE, phi = FALSE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ occasion,
                                                    p = ~ survey.eff + mei),
         group.pars = list(b = FALSE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ occasion,
                                                    p = ~ survey.eff + mei),
         group.pars = list(b = TRUE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ occasion,
                                                    p = ~ survey.eff + mei),
         group.pars = list(b = TRUE, phi = TRUE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ occasion,
                                                    p = ~ survey.eff + mei),
         group.pars = list(b = TRUE, phi = FALSE, p = FALSE),
         df = covs),
    
    ## Models with everything varying by occasion
    ## except for p (by survey effort + lag.mei)
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ occasion,
                                                    p = ~ survey.eff + lag.mei),
         group.pars = list(b = FALSE, phi = FALSE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ occasion,
                                                    p = ~ survey.eff + lag.mei),
         group.pars = list(b = FALSE, phi = FALSE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ occasion,
                                                    p = ~ survey.eff + lag.mei),
         group.pars = list(b = FALSE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ occasion,
                                                    p = ~ survey.eff + lag.mei),
         group.pars = list(b = TRUE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ occasion,
                                                    p = ~ survey.eff + lag.mei),
         group.pars = list(b = TRUE, phi = TRUE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ occasion,
                                                    p = ~ survey.eff + lag.mei),
         group.pars = list(b = TRUE, phi = FALSE, p = FALSE),
         df = covs),
    
    ## Models with everything varying by occasion
    ## except for b (by mei)
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ occasion,
                                                    p = ~ occasion),
         group.pars = list(b = FALSE, phi = FALSE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ occasion,
                                                    p = ~ occasion),
         group.pars = list(b = FALSE, phi = FALSE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ occasion,
                                                    p = ~ occasion),
         group.pars = list(b = FALSE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ occasion,
                                                    p = ~ occasion),
         group.pars = list(b = TRUE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ occasion,
                                                    p = ~ occasion),
         group.pars = list(b = TRUE, phi = TRUE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ occasion,
                                                    p = ~ occasion),
         group.pars = list(b = TRUE, phi = FALSE, p = FALSE),
         df = covs),
    
    ## Models with b varying by mei, phi by occasion, and 
    ## p by survey effort
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ occasion,
                                                    p = ~ survey.eff),
         group.pars = list(b = FALSE, phi = FALSE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ occasion,
                                                    p = ~ survey.eff),
         group.pars = list(b = FALSE, phi = FALSE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ occasion,
                                                    p = ~ survey.eff),
         group.pars = list(b = FALSE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ occasion,
                                                    p = ~ survey.eff),
         group.pars = list(b = TRUE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ occasion,
                                                    p = ~ survey.eff),
         group.pars = list(b = TRUE, phi = TRUE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ occasion,
                                                    p = ~ survey.eff),
         group.pars = list(b = TRUE, phi = FALSE, p = FALSE),
         df = covs),
    
    ## Models with b and p varying by mei, phi by occasion
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ occasion,
                                                    p = ~ mei),
         group.pars = list(b = FALSE, phi = FALSE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ occasion,
                                                    p = ~ mei),
         group.pars = list(b = FALSE, phi = FALSE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ occasion,
                                                    p = ~ mei),
         group.pars = list(b = FALSE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ occasion,
                                                    p = ~ mei),
         group.pars = list(b = TRUE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ occasion,
                                                    p = ~ mei),
         group.pars = list(b = TRUE, phi = TRUE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ occasion,
                                                    p = ~ mei),
         group.pars = list(b = TRUE, phi = FALSE, p = FALSE),
         df = covs),
    
    ## Models with b varying by mei, phi by occasion, and 
    ## p by lag.mei
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ occasion,
                                                    p = ~ lag.mei),
         group.pars = list(b = FALSE, phi = FALSE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ occasion,
                                                    p = ~ lag.mei),
         group.pars = list(b = FALSE, phi = FALSE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ occasion,
                                                    p = ~ lag.mei),
         group.pars = list(b = FALSE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ occasion,
                                                    p = ~ lag.mei),
         group.pars = list(b = TRUE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ occasion,
                                                    p = ~ lag.mei),
         group.pars = list(b = TRUE, phi = TRUE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ occasion,
                                                    p = ~ lag.mei),
         group.pars = list(b = TRUE, phi = FALSE, p = FALSE),
         df = covs),
    
    ## Models with b varying by mei, phi by occasion, and 
    ## p by survey effort + mei
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ occasion,
                                                    p = ~ survey.eff + mei),
         group.pars = list(b = FALSE, phi = FALSE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ occasion,
                                                    p = ~ survey.eff + mei),
         group.pars = list(b = FALSE, phi = FALSE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ occasion,
                                                    p = ~ survey.eff + mei),
         group.pars = list(b = FALSE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ occasion,
                                                    p = ~ survey.eff + mei),
         group.pars = list(b = TRUE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ occasion,
                                                    p = ~ survey.eff + mei),
         group.pars = list(b = TRUE, phi = TRUE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ occasion,
                                                    p = ~ survey.eff + mei),
         group.pars = list(b = TRUE, phi = FALSE, p = FALSE),
         df = covs),
    
    ## Models with b varying by mei, phi by occasion, and 
    ## p by survey effort + lag.mei
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ occasion,
                                                    p = ~ survey.eff + lag.mei),
         group.pars = list(b = FALSE, phi = FALSE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ occasion,
                                                    p = ~ survey.eff + lag.mei),
         group.pars = list(b = FALSE, phi = FALSE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ occasion,
                                                    p = ~ survey.eff + lag.mei),
         group.pars = list(b = FALSE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ occasion,
                                                    p = ~ survey.eff + lag.mei),
         group.pars = list(b = TRUE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ occasion,
                                                    p = ~ survey.eff + lag.mei),
         group.pars = list(b = TRUE, phi = TRUE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ occasion,
                                                    p = ~ survey.eff + lag.mei),
         group.pars = list(b = TRUE, phi = FALSE, p = FALSE),
         df = covs),
    
    ## Models with b varying by lag.mei, phi and p by occasion
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ occasion,
                                                    p = ~ occasion),
         group.pars = list(b = FALSE, phi = FALSE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ occasion,
                                                    p = ~ occasion),
         group.pars = list(b = FALSE, phi = FALSE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ occasion,
                                                    p = ~ occasion),
         group.pars = list(b = FALSE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ occasion,
                                                    p = ~ occasion),
         group.pars = list(b = TRUE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ occasion,
                                                    p = ~ occasion),
         group.pars = list(b = TRUE, phi = TRUE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ occasion,
                                                    p = ~ occasion),
         group.pars = list(b = TRUE, phi = FALSE, p = FALSE),
         df = covs),
    
    ## Models with b varying by lag.mei, phi by occasion
    ## p by survey effort
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ occasion,
                                                    p = ~ survey.eff),
         group.pars = list(b = FALSE, phi = FALSE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ occasion,
                                                    p = ~ survey.eff),
         group.pars = list(b = FALSE, phi = FALSE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ occasion,
                                                    p = ~ survey.eff),
         group.pars = list(b = FALSE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ occasion,
                                                    p = ~ survey.eff),
         group.pars = list(b = TRUE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ occasion,
                                                    p = ~ survey.eff),
         group.pars = list(b = TRUE, phi = TRUE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ occasion,
                                                    p = ~ survey.eff),
         group.pars = list(b = TRUE, phi = FALSE, p = FALSE),
         df = covs),

    ## Models with b varying by lag.mei, phi by occasion
    ## p by mei
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ occasion,
                                                    p = ~ mei),
         group.pars = list(b = FALSE, phi = FALSE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ occasion,
                                                    p = ~ mei),
         group.pars = list(b = FALSE, phi = FALSE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ occasion,
                                                    p = ~ mei),
         group.pars = list(b = FALSE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ occasion,
                                                    p = ~ mei),
         group.pars = list(b = TRUE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ occasion,
                                                    p = ~ mei),
         group.pars = list(b = TRUE, phi = TRUE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ occasion,
                                                    p = ~ mei),
         group.pars = list(b = TRUE, phi = FALSE, p = FALSE),
         df = covs),
    
    ## Models with b and p varying by lag.mei, phi by occasion
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ occasion,
                                                    p = ~ lag.mei),
         group.pars = list(b = FALSE, phi = FALSE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ occasion,
                                                    p = ~ lag.mei),
         group.pars = list(b = FALSE, phi = FALSE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ occasion,
                                                    p = ~ lag.mei),
         group.pars = list(b = FALSE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ occasion,
                                                    p = ~ lag.mei),
         group.pars = list(b = TRUE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ occasion,
                                                    p = ~ lag.mei),
         group.pars = list(b = TRUE, phi = TRUE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ occasion,
                                                    p = ~ lag.mei),
         group.pars = list(b = TRUE, phi = FALSE, p = FALSE),
         df = covs),
    
    ## Models with b varying by lag.mei, phi by occasion
    ## p by survey effort + mei
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ occasion,
                                                    p = ~ survey.eff + mei),
         group.pars = list(b = FALSE, phi = FALSE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ occasion,
                                                    p = ~ survey.eff + mei),
         group.pars = list(b = FALSE, phi = FALSE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ occasion,
                                                    p = ~ survey.eff + mei),
         group.pars = list(b = FALSE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ occasion,
                                                    p = ~ survey.eff + mei),
         group.pars = list(b = TRUE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ occasion,
                                                    p = ~ survey.eff + mei),
         group.pars = list(b = TRUE, phi = TRUE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ occasion,
                                                    p = ~ survey.eff + mei),
         group.pars = list(b = TRUE, phi = FALSE, p = FALSE),
         df = covs),
    
    ## Models with b varying by lag.mei, phi by occasion
    ## p by survey effort + lag.mei
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ occasion,
                                                    p = ~ survey.eff + lag.mei),
         group.pars = list(b = FALSE, phi = FALSE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ occasion,
                                                    p = ~ survey.eff + lag.mei),
         group.pars = list(b = FALSE, phi = FALSE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ occasion,
                                                    p = ~ survey.eff + lag.mei),
         group.pars = list(b = FALSE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ occasion,
                                                    p = ~ survey.eff + lag.mei),
         group.pars = list(b = TRUE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ occasion,
                                                    p = ~ survey.eff + lag.mei),
         group.pars = list(b = TRUE, phi = TRUE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ occasion,
                                                    p = ~ survey.eff + lag.mei),
         group.pars = list(b = TRUE, phi = FALSE, p = FALSE),
         df = covs),
    
    ## Models with b varying by occasion, constant phi
    ## p by occasion
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ 1,
                                                    p = ~ occasion),
         group.pars = list(b = FALSE, phi = FALSE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ 1,
                                                    p = ~ occasion),
         group.pars = list(b = FALSE, phi = FALSE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ 1,
                                                    p = ~ occasion),
         group.pars = list(b = FALSE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ 1,
                                                    p = ~ occasion),
         group.pars = list(b = TRUE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ 1,
                                                    p = ~ occasion),
         group.pars = list(b = TRUE, phi = TRUE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ 1,
                                                    p = ~ occasion),
         group.pars = list(b = TRUE, phi = FALSE, p = FALSE),
         df = covs),
    
    ## Models with b varying by occasion, constant phi
    ## p by survey effort
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ 1,
                                                    p = ~ survey.eff),
         group.pars = list(b = FALSE, phi = FALSE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ 1,
                                                    p = ~ survey.eff),
         group.pars = list(b = FALSE, phi = FALSE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ 1,
                                                    p = ~ survey.eff),
         group.pars = list(b = FALSE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ 1,
                                                    p = ~ survey.eff),
         group.pars = list(b = TRUE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ 1,
                                                    p = ~ survey.eff),
         group.pars = list(b = TRUE, phi = TRUE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ 1,
                                                    p = ~ survey.eff),
         group.pars = list(b = TRUE, phi = FALSE, p = FALSE),
         df = covs),

    ## Models with b varying by occasion, constant phi
    ## p by mei
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ 1,
                                                    p = ~ mei),
         group.pars = list(b = FALSE, phi = FALSE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ 1,
                                                    p = ~ mei),
         group.pars = list(b = FALSE, phi = FALSE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ 1,
                                                    p = ~ mei),
         group.pars = list(b = FALSE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ 1,
                                                    p = ~ mei),
         group.pars = list(b = TRUE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ 1,
                                                    p = ~ mei),
         group.pars = list(b = TRUE, phi = TRUE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ 1,
                                                    p = ~ mei),
         group.pars = list(b = TRUE, phi = FALSE, p = FALSE),
         df = covs),
    
    ## Models with b varying by occasion, constant phi
    ## p by lag.mei
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ 1,
                                                    p = ~ lag.mei),
         group.pars = list(b = FALSE, phi = FALSE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ 1,
                                                    p = ~ lag.mei),
         group.pars = list(b = FALSE, phi = FALSE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ 1,
                                                    p = ~ lag.mei),
         group.pars = list(b = FALSE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ 1,
                                                    p = ~ lag.mei),
         group.pars = list(b = TRUE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ 1,
                                                    p = ~ lag.mei),
         group.pars = list(b = TRUE, phi = TRUE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ 1,
                                                    p = ~ lag.mei),
         group.pars = list(b = TRUE, phi = FALSE, p = FALSE),
         df = covs),
    
    ## Models with b varying by occasion, constant phi
    ## p by survey effort + mei
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ 1,
                                                    p = ~ survey.eff + mei),
         group.pars = list(b = FALSE, phi = FALSE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ 1,
                                                    p = ~ survey.eff + mei),
         group.pars = list(b = FALSE, phi = FALSE, p = TRUE),
         df = covs),

    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ 1,
                                                    p = ~ survey.eff + mei),
         group.pars = list(b = FALSE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ 1,
                                                    p = ~ survey.eff + mei),
         group.pars = list(b = TRUE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ 1,
                                                    p = ~ survey.eff + mei),
         group.pars = list(b = TRUE, phi = TRUE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ 1,
                                                    p = ~ survey.eff + mei),
         group.pars = list(b = TRUE, phi = FALSE, p = FALSE),
         df = covs),
    
    ## Models with b varying by occasion, constant phi
    ## p by survey effort + lag.mei
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ 1,
                                                    p = ~ survey.eff + lag.mei),
         group.pars = list(b = FALSE, phi = FALSE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ 1,
                                                    p = ~ survey.eff + lag.mei),
         group.pars = list(b = FALSE, phi = FALSE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ 1,
                                                    p = ~ survey.eff + lag.mei),
         group.pars = list(b = FALSE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ 1,
                                                    p = ~ survey.eff + lag.mei),
         group.pars = list(b = TRUE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ 1,
                                                    p = ~ survey.eff + lag.mei),
         group.pars = list(b = TRUE, phi = TRUE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ occasion,
                                                    phi = ~ 1,
                                                    p = ~ survey.eff + lag.mei),
         group.pars = list(b = TRUE, phi = FALSE, p = FALSE),
         df = covs),
    
    ## Models with b varying by mei, constant phi
    ## p by occasion
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ 1,
                                                    p = ~ occasion),
         group.pars = list(b = FALSE, phi = FALSE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ 1,
                                                    p = ~ occasion),
         group.pars = list(b = FALSE, phi = FALSE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ 1,
                                                    p = ~ occasion),
         group.pars = list(b = FALSE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ 1,
                                                    p = ~ occasion),
         group.pars = list(b = TRUE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ 1,
                                                    p = ~ occasion),
         group.pars = list(b = TRUE, phi = TRUE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ 1,
                                                    p = ~ occasion),
         group.pars = list(b = TRUE, phi = FALSE, p = FALSE),
         df = covs),
    
    ## Models with b varying by mei, constant phi
    ## p by survey effort
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ 1,
                                                    p = ~ survey.eff),
         group.pars = list(b = FALSE, phi = FALSE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ 1,
                                                    p = ~ survey.eff),
         group.pars = list(b = FALSE, phi = FALSE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ 1,
                                                    p = ~ survey.eff),
         group.pars = list(b = FALSE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ 1,
                                                    p = ~ survey.eff),
         group.pars = list(b = TRUE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ 1,
                                                    p = ~ survey.eff),
         group.pars = list(b = TRUE, phi = TRUE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ 1,
                                                    p = ~ survey.eff),
         group.pars = list(b = TRUE, phi = FALSE, p = FALSE),
         df = covs),
    
    ## Models with b and p varying by mei, constant phi
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ 1,
                                                    p = ~ mei),
         group.pars = list(b = FALSE, phi = FALSE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ 1,
                                                    p = ~ mei),
         group.pars = list(b = FALSE, phi = FALSE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ 1,
                                                    p = ~ mei),
         group.pars = list(b = FALSE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ 1,
                                                    p = ~ mei),
         group.pars = list(b = TRUE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ 1,
                                                    p = ~ mei),
         group.pars = list(b = TRUE, phi = TRUE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ 1,
                                                    p = ~ mei),
         group.pars = list(b = TRUE, phi = FALSE, p = FALSE),
         df = covs),
    
    ## Models with b varying by mei, constant phi
    ## p by lag.mei
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ 1,
                                                    p = ~ lag.mei),
         group.pars = list(b = FALSE, phi = FALSE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ 1,
                                                    p = ~ lag.mei),
         group.pars = list(b = FALSE, phi = FALSE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ 1,
                                                    p = ~ lag.mei),
         group.pars = list(b = FALSE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ 1,
                                                    p = ~ lag.mei),
         group.pars = list(b = TRUE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ 1,
                                                    p = ~ lag.mei),
         group.pars = list(b = TRUE, phi = TRUE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ 1,
                                                    p = ~ lag.mei),
         group.pars = list(b = TRUE, phi = FALSE, p = FALSE),
         df = covs),
    
  
    ## Models with b varying by mei, constant phi
    ## p by survey effort + mei
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ 1,
                                                    p = ~ survey.eff + mei),
         group.pars = list(b = FALSE, phi = FALSE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ 1,
                                                    p = ~ survey.eff + mei),
         group.pars = list(b = FALSE, phi = FALSE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ 1,
                                                    p = ~ survey.eff + mei),
         group.pars = list(b = FALSE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ 1,
                                                    p = ~ survey.eff + mei),
         group.pars = list(b = TRUE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ 1,
                                                    p = ~ survey.eff + mei),
         group.pars = list(b = TRUE, phi = TRUE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ 1,
                                                    p = ~ survey.eff + mei),
         group.pars = list(b = TRUE, phi = FALSE, p = FALSE),
         df = covs),
    
    ## Models with b varying by mei, constant phi
    ## p by survey effort + lag.mei
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ 1,
                                                    p = ~ survey.eff + lag.mei),
         group.pars = list(b = FALSE, phi = FALSE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ 1,
                                                    p = ~ survey.eff + lag.mei),
         group.pars = list(b = FALSE, phi = FALSE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ 1,
                                                    p = ~ survey.eff + lag.mei),
         group.pars = list(b = FALSE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ 1,
                                                    p = ~ survey.eff + lag.mei),
         group.pars = list(b = TRUE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ 1,
                                                    p = ~ survey.eff + lag.mei),
         group.pars = list(b = TRUE, phi = TRUE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ mei,
                                                    phi = ~ 1,
                                                    p = ~ survey.eff + lag.mei),
         group.pars = list(b = TRUE, phi = FALSE, p = FALSE),
         df = covs),
    
    ## Models with b varying by lag.mei, constant phi
    ## p by occasion
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ 1,
                                                    p = ~ occasion),
         group.pars = list(b = FALSE, phi = FALSE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ 1,
                                                    p = ~ occasion),
         group.pars = list(b = FALSE, phi = FALSE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ 1,
                                                    p = ~ occasion),
         group.pars = list(b = FALSE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ 1,
                                                    p = ~ occasion),
         group.pars = list(b = TRUE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ 1,
                                                    p = ~ occasion),
         group.pars = list(b = TRUE, phi = TRUE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ 1,
                                                    p = ~ occasion),
         group.pars = list(b = TRUE, phi = FALSE, p = FALSE),
         df = covs),
    
    ## Models with b varying by lag.mei, constant phi
    ## p by survey effort
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ 1,
                                                    p = ~ survey.eff),
         group.pars = list(b = FALSE, phi = FALSE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ 1,
                                                    p = ~ survey.eff),
         group.pars = list(b = FALSE, phi = FALSE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ 1,
                                                    p = ~ survey.eff),
         group.pars = list(b = FALSE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ 1,
                                                    p = ~ survey.eff),
         group.pars = list(b = TRUE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ 1,
                                                    p = ~ survey.eff),
         group.pars = list(b = TRUE, phi = TRUE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ 1,
                                                    p = ~ survey.eff),
         group.pars = list(b = TRUE, phi = FALSE, p = FALSE),
         df = covs),
    
    ## Models with b varying by lag.mei, constant phi
    ## p by mei
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ 1,
                                                    p = ~ mei),
         group.pars = list(b = FALSE, phi = FALSE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ 1,
                                                    p = ~ mei),
         group.pars = list(b = FALSE, phi = FALSE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ 1,
                                                    p = ~ mei),
         group.pars = list(b = FALSE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ 1,
                                                    p = ~ mei),
         group.pars = list(b = TRUE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ 1,
                                                    p = ~ mei),
         group.pars = list(b = TRUE, phi = TRUE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ 1,
                                                    p = ~ mei),
         group.pars = list(b = TRUE, phi = FALSE, p = FALSE),
         df = covs),
    
    ## Models with b varying by lag.mei, constant phi
    ## p by lag.mei
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ 1,
                                                    p = ~ lag.mei),
         group.pars = list(b = FALSE, phi = FALSE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ 1,
                                                    p = ~ lag.mei),
         group.pars = list(b = FALSE, phi = FALSE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ 1,
                                                    p = ~ lag.mei),
         group.pars = list(b = FALSE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ 1,
                                                    p = ~ lag.mei),
         group.pars = list(b = TRUE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ 1,
                                                    p = ~ lag.mei),
         group.pars = list(b = TRUE, phi = TRUE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ 1,
                                                    p = ~ lag.mei),
         group.pars = list(b = TRUE, phi = FALSE, p = FALSE),
         df = covs),
    
    ## Models with b varying by lag.mei, constant phi
    ## p by survey effort + mei
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ 1,
                                                    p = ~ survey.eff + mei),
         group.pars = list(b = FALSE, phi = FALSE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ 1,
                                                    p = ~ survey.eff + mei),
         group.pars = list(b = FALSE, phi = FALSE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ 1,
                                                    p = ~ survey.eff + mei),
         group.pars = list(b = FALSE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ 1,
                                                    p = ~ survey.eff + mei),
         group.pars = list(b = TRUE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ 1,
                                                    p = ~ survey.eff + mei),
         group.pars = list(b = TRUE, phi = TRUE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ 1,
                                                    p = ~ survey.eff + mei),
         group.pars = list(b = TRUE, phi = FALSE, p = FALSE),
         df = covs),
    
    ## Models with b varying by lag.mei, constant phi
    ## p by survey effort + lag.mei
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ 1,
                                                    p = ~ survey.eff + lag.mei),
         group.pars = list(b = FALSE, phi = FALSE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ 1,
                                                    p = ~ survey.eff + lag.mei),
         group.pars = list(b = FALSE, phi = FALSE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ 1,
                                                    p = ~ survey.eff + lag.mei),
         group.pars = list(b = FALSE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ 1,
                                                    p = ~ survey.eff + lag.mei),
         group.pars = list(b = TRUE, phi = TRUE, p = TRUE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ 1,
                                                    p = ~ survey.eff + lag.mei),
         group.pars = list(b = TRUE, phi = TRUE, p = FALSE),
         df = covs),
    list(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                    phi = ~ 1,
                                                    p = ~ survey.eff + lag.mei),
         group.pars = list(b = TRUE, phi = FALSE, p = FALSE),
         df = covs)
)

## fits all possible models on the list above ####
fits <- par.fit.popan(3, arg.list = args)
save(fits, file = "2020.10.15_fits.RData")
load("2020.10.15_fits.RData")


# list AICs of all possible models
sapply(fits, AIC)
# sort AICs of all models from the lowest to the highest
order(sapply(fits, AIC))

do.call("fit.popan", args[[186]])

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
## Summary for estimated parameters.
summary.popan(boot.best.fit, function(fit) fit$fit$par)
## Summary for ENs for females.
summary.popan(boot.best.fit, function(fit) fit$ENs[, 1])
## Summary for ENs for males.
summary.popan(boot.best.fit, function(fit) fit$ENs[, 2])
## Summary for ENs for both groups.
summary.popan(boot.best.fit, function(fit) fit$ENs[, 1] + fit$ENs[, 2])
## Summary comparing first and last years (1 = females, and 2= males)
summary.popan(boot.best.fit, function(fit) fit$ENs[11, ] - fit$ENs[1, ])
## Summary comparing females and males
summary.popan(boot.best.fit, function(fit) fit$ENs[, 1] - fit$ENs[, 2])

## For p-value comparing first and last year for females:
## (1) calculate t-test statistic: estimate/std.error from the summary
372.6852/34.53555
## (2) compare to a standard normal distribution to get the p-value
## don't forget the negative value before 'estimate/std.error
## and multiply by 2 for a two-tailed test:
2*pnorm(-372.6852/34.53555)

## For p-value comparing first and last year for males:
## (1) calculate t-test statistic: estimate/std.error:
143.5522/63.11903
## (2) compare to a standard normal distribution to get the p-value
## and multiply by 2 for a two-tailed test:
2*pnorm(-143.5522/63.11903)

## For p-value comparing females and males in 2019 
## (if females were more abundant than males)
## (1) calculate t-test statistic: estimate/standard error
240.26506/59.06843
## (2) compare to a standard normal distribution to get the p-value
## and multiply by 2 for a two-tailed test:
2*pnorm(-240.26506/59.06843)

## Comparing males to females in terms of survival.
summary.popan(boot.best.fit, function(fit) fit$phis[1:2, 1] - fit$phis[1:2, 2])
# 2*pnorm(-abs(Estimate/Std Error)): use Estimate and Std Error from the summary above
# females [1,] and males [2,]
2*pnorm(-abs(0.05945483/0.02612161))
## Restricted fit with male and female survival the same.
phi.restricted.fit <- fit.popan(caplist = captlist[1:2], model.list = list(b = ~ lag.mei,
                                                                       phi = ~ 1,
                                                                       p = ~ occasion),
                            group.pars = list(b = TRUE, phi = TRUE, p = FALSE),
                            df = covs)
## Likelihood ratio test statistics is twice difference in log-likelihoods.
phi.lrts <- 2*(phi.restricted.fit$fit$objective - best.fit$fit$objective)
npar.diff <- length(best.fit$fit$par) - length(phi.restricted.fit$fit$par)
## P-value is calculated by comparing LRTS to a chi-squared distribution with npar.diff degrees of freedom.
1 - pchisq(phi.lrts, npar.diff)

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
1 - pchisq(b.lrts, b.npar.diff) ## result: 0


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
save(boot.best.fit, file = "2020.10.15_boot-best-fit-manta.RData")


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
plot(covs$unlag.mei)

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


