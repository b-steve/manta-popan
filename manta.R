## Sourcing model fitting functions.
source("Rfunc.R")
## Reading in data.
load("manta.RData")



captlist.sim <- simGroups.popanGeneral.func(k=11,
                                            birthfunc = immigrationElNino.func,
                                            Nsvec=c(1000, 1200),
                                            ## Different survival by group and over time: each phivec should have k-1=10 elts:
                                            phiList=list(c(rep(0.8, 6), 0.4, 0.4, rep(0.8, 2)), c(rep(0.9, 6), 0.5, 0.5, rep(0.9, 2))),
                                            birthparList=list(bpar1 = c(0.18, 0.1)),
                                            ## just one list element means these parameters will be used for both groups
                                            pList= list(rep(0.3, 11), rep(0.2, 11)))  ## two list elements => different by group


birth.elnino.mei <- function(bpars = c(b1 = 0, b2 = 0)){
    b1 <- bpars[1]
    b2 <- bpars[2]
    ## Note that we do not need rho for the final year.
    eta <- b1 + b2*covs$mei[-length(covs$mei)]
    exp(eta)
}

## A popan-General model:
popanGeneral.fit.func(captlist[1:2], k=11, birthfunc = birth.elnino.mei,
                      model=list(
                          gp1=c("Ns.1",
                                "b1.1","b2.1",
                                c(rep("phi1.1", 6), "phi7.1","phi7.1", rep("phi1.1", 2)),
                                rep("p1.1", 11)),
                          gp1=c("Ns.2",
                                "b1.1","b2.1",
                                c(rep("phi1.2", 6), "phi7.2","phi7.2", rep("phi1.2", 2)),
                                rep("p1.2", 11))
                      ),
                      ## start values: this named vector can be entered in any order
                      startvec=c(Ns.1=1000, Ns.2=800, b1.1= log(0.18), b2.1=0, phi1.1=0.9, phi7.1=0.2,
                                 phi1.2=0.6, phi7.2=0.4, p1.1=0.1, p1.2=0.3), printit = TRUE)

