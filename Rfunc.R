if(F){
    ## Code originally written by Robin Aldridge-Sutton, adjusted by Rachel Fewster to facilitate fitting
    ## for Manta rays with multiple groups (e.g. sex) and covariates.

    ## -------------------------------
    ## ORGANISATION:
    ## -------------------------------
    ## I've separated the coding for POPAN-lambda models from that for POPAN-general models.
    ## This is just for illustrative purposes so that you can see what's going on in the POPAN-lambda model
    ## without the extra skullduggery of the general formulation.  Because POPAN-lambda is a special case
    ## of POPAN-general, the lambda-specific functions can be discarded ultimately: there's no need to
    ## maintain separate copies of the code in different functions.
    ##
    ## Inventory:
    ## POPAN-lambda simulation and fitting functions: simGroups.popanLambda.func, popanLambda.fit.func
    ## POPAN-general simulation and fitting functions: simGroups.popanGeneral.func, popanGeneral.fit.func
    ## Wrapper functions to simulate and fit many times: popanLambda.wrap, popanGeneral.wrap
    ## Functions for calculating pentvec under the two model types: pentLambda.func, pentGeneral.func
    ## Robin's functions for the likelihood of any model type: popansummary.func, chi.func, psi.func, ptheta.func
    ## Example function to show how to incorporate covariates into birth/immigration:  immigrationElNino.func

    ## -------------------------------
    ## WHAT'S MISSING:
    ## -------------------------------
    ## In the interests of code readability, I haven't included various things that usually would be included,
    ## such as finding start-vectors for optimization, estimating variance, plotting results, and error-catching
    ## during repeated fits. I figure Ben will likely have his preferred ways of doing all the above - which might
    ## include translating into TMB, resolving the start-vector and variance estimation in one go.
    ## The aim of this code is to establish the population dynamics modelling and create a general structure
    ## for incorporating covariates. However there is plenty more code for pretty-much everything if you want it.

    ## --------------------------------
    ## POPAN-lambda DEMO:
    ## --------------------------------
    ## Source the code:
    source("Rfunc.R")

    ## Simulate two groups with different superpopulation sizes, and capture
    ## probabilities, but same survival probability (phi) and per-capita birth rate (rho):
    ## note that this model has been reparametrized in terms of (phi, rho) rather than (phi, lambda):
    ## see comments below in simGroups.popanLambda.func.  It is still a "POPAN-lambda" model because
    ## phi and rho remain constant through time. The relationship is lambda = phi + rho.
    caplist <- simGroups.popanLambda.func(k=11,
                                          Nsvec=c(1000, 1200),
                                          phiList=0.85,
                                          rhoList=0.18,
                                          pList= list(rep(0.3, 11), rep(0.2, 11)))

    ## Fit popan-Lambda model with single capture-probability for each group, common rho and phi between
    ## groups. This could be suitable if the two groups are (male, female) and abundance and capture probability
    ## differ for the two sexes.  Observe how we specify the model, with rho and phi common between the two
    ## groups, and constant capture probabilities within group: pt.1=p1.1 and pt.2=p1.2 for all t.
    ## Consistently, anything ending in ".1" means group 1 and ".2" means group 2.
    popanLambda.fit.func(caplist, k=11, model=list(
                                            gp1=c("Ns.1", "rho.1", "phi.1", rep("p1.1", 11)),
                                            gp2=c("Ns.2", "rho.1", "phi.1", rep("p1.2", 11))),
                         ## start from the generating values: the optimization needs good start-points.
                         startvec=c(Ns.1=1000, Ns.2=1200, rho.1=0.18, phi.1=0.85, p1.1=0.3, p1.2=0.2), printit = TRUE)

    ## Repeat the above 200 times and plot estimates:
    sim.res <- popanLambda.wrap(Nsim=200)
    ## Plot Ns estimates:
    boxplot(sim.res[grep("Ns", names(sim.res))], col="lightgoldenrod")
    abline(h=c(1000, 1200), col="blue")  ## Check against true values
    ## Plot all other estimates:
    boxplot(sim.res[-grep("Ns", names(sim.res))], col="orchid")
    abline(h=c(0.18, 0.85, 0.3, 0.2), col="grey")  ## Check against true values

    ## ---------------------------------
    ## POPAN-General DEMO:
    ## ---------------------------------
    ## For example, use the simple function immigrationElNino.func to establish higher per-capita "birth" rates in
    ## two years, specifically 2014 and 2015, which respectively apply to immigrants arriving in 2015 and 2016.
    ## The code is set up so that any function can be used here, as long as it takes a vector of parameters
    ## as inputs and returns a vector of overall per-capita birthrates of length k-1.

    ## Note that POPAN-General allows for different survival probabilities by year, as well as allowing
    ## birth/immigration to be established by a generic covariate-based function.

    ## First simulate two groups with different superpopulation sizes, and capture probabilities, but the same
    ## (and constant) survival probability phi, and the same per-capita parameters to the birth-function (b):
    ## the birth-function is a supplied argument.
    caplist <- simGroups.popanGeneral.func(k=11,
                                           birthfunc = immigrationElNino.func,
                                           Nsvec=c(1000, 1200),
                                           phiList=0.85,  ## specifying one phi means it will be cycled for all times and groups
                                           birthparList=list(bpar1 = c(0.18, 0.1)),
                                           ## just one list element means these parameters will be used for both groups
                                           pList= list(rep(0.3, 11), rep(0.2, 11)))  ## two list elements => different by group

    ## Fit popan-General model with single capture-probability for each group, and common birth parameters
    ## and phi between groups.  Parameters must be entered in the right order in the model argument below,
    ## and all parameters must be accounted for: Ns, bpars (as many as needed), k-1 phipars, k ppars.
    ## Consistently, anything ending in ".1" means group 1 and ".2" means group 2.
    popanGeneral.fit.func(caplist, k=11, birthfunc = immigrationElNino.func, model=list(
                                                                                 gp1=c("Ns.1", "b1.1", "b2.1", rep("phi1.1", 10), rep("p1.1", 11)),
                                                                                 gp2=c("Ns.2", "b1.1", "b2.1", rep("phi1.1", 10), rep("p1.2", 11))),
                          ## start from the generating values: the optimization needs good start-points.
                          startvec=c(Ns.1=1000, Ns.2=1200, b1.1=0.18, b2.1=0.1, phi1.1=0.85,
                                     p1.1=0.3, p1.2=0.2), printit = TRUE)

    ## Repeat the above 200 times and plot estimates:
    sim.res <- popanGeneral.wrap(Nsim=200)
    ## Plot Ns estimates:
    boxplot(sim.res[grep("Ns", names(sim.res))], col="lightgoldenrod")
    abline(h=c(1000, 1200), col="blue")  ## Check against true values
    ## Plot all other estimates:
    boxplot(sim.res[-grep("Ns", names(sim.res))], col="orchid")
    abline(h=c(0.18, 0.1, 0.85, 0.3, 0.2), col="grey")  ## Check against true values

    ## -------------------------------------------------------------------------------------------------------------------
    ## Demonstrate survival varying by both time and group, as well as the ElNino immigration model:
    caplist <- simGroups.popanGeneral.func(k=11,
                                           birthfunc = immigrationElNino.func,
                                           Nsvec=c(1000, 1200),
                                           ## Different survival by group and over time: each phivec should have k-1=10 elts:
                                           phiList=list(c(rep(0.8, 6), 0.4, 0.4, rep(0.8, 2)), c(rep(0.9, 6), 0.5, 0.5, rep(0.9, 2))),
                                           birthparList=list(bpar1 = c(0.18, 0.1)),
                                           ## just one list element means these parameters will be used for both groups
                                           pList= list(rep(0.3, 11), rep(0.2, 11)))  ## two list elements => different by group

    ## Fit the corresponding popan-General model:
    popanGeneral.fit.func(caplist, k=11, birthfunc = immigrationElNino.func,
                          model=list(
                              gp1=c("Ns.1","b1.1","b2.1",c(rep("phi1.1", 6), "phi7.1","phi7.1", rep("phi1.1", 2)), rep("p1.1", 11)),
                              gp1=c("Ns.2","b1.1","b2.1",c(rep("phi1.2", 6), "phi7.2","phi7.2", rep("phi1.2", 2)), rep("p1.2", 11))
                          ),
                          ## start from the generating values: this named vector can be entered in any order
                          startvec=c(Ns.1=1000, Ns.2=1200, b1.1=0.18, b2.1=0.1, phi1.1=0.8, phi7.1=0.4,
                                     phi1.2=0.9, phi7.2=0.5, p1.1=0.3, p1.2=0.2), printit = TRUE)



}

###############################################################
## FUNCTIONS BELOW
###############################################################

pentLambda.func <- function(rho, phi, k){
    ## pentLambda.func is the p-ent function for the popan-Lambda model, parametrized via rho=lambda-phi:
    ## This function assumes there are surveys every year (no survey gaps).
    pentvec <- c(1, rho * (rho+phi)^(0:(k - 2)))
    pentvec / sum(pentvec)
}


###############################################################

pentGeneral.func <- function(rhovec, phivec, k){
    ## pentGeneral.func takes a vector of birth-rates, rhovec, and survival probabilities, phivec,
    ## and calculates the corresponding pentvec.
    ## Each of rhovec and phivec have length k-1.
    ## rhovec[t] is the apparent birth-rate from time t to time t+1, including both births and immigration;
    ## so the expected number of animals joining the population for the first time at time t+1 is
    ## rhovec[t]*Nt.
    ## phivec[t] is the apparent survival from time t to time t+1, allowing for both deaths and emigration.
    ## So the expected number in the population at time t+1 is:
    ## E(N_{t+1} | Nt) = phi[t]*Nt + rho[t]*Nt = (phi[t]+rho[t])*Nt.
    ## Averaging over Nt gives EN(t+1) = (phi[t] + rho[t])*ENt.  (1)
    ##
    ## Now EN(t+1) is also expressed in terms of the superpopulation size, Ns, and pent[t+1], as follows:
    ## EN(t+1) = phi[t]*ENt + pent[t+1]*Ns. (2)
    ##
    ## Equating (1) and (2) gives rho[t]*ENt = pent[t+1]*Ns, so pent[t+1] = rho[t]*ENt / Ns.
    ##
    ## The initialisation is pent1 = EN1/Ns, then:
    ## pent2=rho1*EN1/Ns = rho1*pent1
    ## pent3=rho2*EN2/Ns = rho2*(phi1+rho1)*EN1/Ns = rho2*(phi1+rho1)*pent1
    ## pent4=rho3*EN3/Ns = rho3*(phi2+rho2)*EN2/Ns = rho3*(phi2+rho2)*(phi1+rho1)*pent1
    ## etc.
    ## In general, pent[t+1] = rho[t]*prod( (phi+rho)[1:t-1] ) * pent1.
    ## So pent[2:k] = rho[1:k-1]*cumprod(c(1, (phi+rho)[1:k-2])) * pent1.
    ##
    ## For the POPAN-lambda model, set phi and rho as constants: then pent[t+1]=rho*(phi+rho)^(t-1)*pent1.
    pentvec <- c(1, rhovec*cumprod(c(1, (phivec+rhovec)[1:(k-2)])))
    pentvec / sum(pentvec)
}


###############################################################
## Function to summarise data for popan model
popansummary.func <- function(det.dat) {
    ## Store data as matrix not data frame
    det.dat <- as.matrix(det.dat)

    ## Find number of capture histories
    nhist <- nrow(det.dat)

    ## Find number of surveys
    k <- ncol(det.dat)

    ## Create empty objects
    non.caps.mat <- matrix(0, nrow = nhist, ncol = k)
    survive.mat <- matrix(0, nrow = nhist, ncol = k)
    first.obs <- numeric(nhist)
    last.obs <- numeric(nhist)

    ## Loop over capture histories in data
    for(elt in 1:nhist){
        ## Find first and last captures
        which.1 <- which(det.dat[elt, ] == 1)
        first <- min(which.1)
        last <- max(which.1)
        first.obs[elt] <- first
        last.obs[elt] <- last

        ## Find non-captures and survival between first and last captures
        if(last > first){
            non.caps.mat[elt, first:last] <- 1 - det.dat[elt, first:last]
            survive.mat[elt, first:(last - 1)] <- 1
        }
    }

    ## Find summaries
    caps <- colSums(det.dat)
    non.caps <- colSums(non.caps.mat)
    first.tab <- table(factor(first.obs, levels=1:k))
    last.tab <- table(factor(last.obs, levels=1:k))
    survives <- colSums(survive.mat)[1:(k-1)]

    ## Combine and return summaries
    summaries <- data.frame(
        first.tab = as.vector(first.tab),
        caps = caps,
        non.caps = non.caps,
        survives = c(survives, NA),
        last.tab = as.vector(last.tab)
    )
    attributes(summaries)$nhist <- nhist
    summaries
}


##############################################################
## Chi function
chi.func <- function(phivec, pvec, k) {
    chivec <- numeric(k)
    chivec[k] <- 1
    for(i in (k - 1):1) {
        chivec[i] <- 1 - phivec[i] + phivec[i] * (1 - pvec[i + 1]) * chivec[i + 1]
    }
    chivec
}

##############################################################
## Psi function
psi.func <- function(pentvec, phivec, pvec, k) {
    psivec <- numeric(k)
    psivec[1] <- pentvec[1]
    for(i in 2:k) {
        psivec[i] <- psivec[i - 1] * (1 - pvec[i - 1]) * phivec[i-1] + pentvec[i]
    }
    psivec
}


##############################################################
## ptheta function
ptheta.func <- function(pentvec, pvec, chivec) {
    1 - sum(pentvec * (1 - pvec) * chivec)
}

###############################################################
simGroups.popanLambda.func <- function(k=11,
                                       Nsvec=c(1000, 1200),
                                       phiList=list(0.85, 0.85),
                                       rhoList=list(0.18, 0.18),
                                       pList= list(rep(0.3, k), rep(0.2, k)))
{
    ## simGroups.popanLambda.func 11/5/2020
    ## Simulate from the POPAN-lambda model for multiple groups (e.g. males and females).
    ## Template by Robin Aldridge-Sutton, recast to loop through multiple groups by Rachel Fewster.
    ##
    ## Note: although it's called the POPAN-lambda model, I've switched the parametrization to rho=lambda-phi,
    ## where rho represents the apparent birth-rate per capita within the group. For example,
    ## if the group population size in year t is Nt, then the number of new animals entering this group by birth
    ## and immigration in year t+1 is rho*Nt.  The net population size in year t+1 is:
    ## E(N_t+1 | Nt) = phi*Nt + rho*Nt = (phi+rho)*Nt = lambda*Nt as usual.
    ## The reason for this switch is that it might be reasonable for the different groups to have the same birth-rate
    ## but different survival rates. This would require 4 parameters in the (phi, lambda) parametrization,
    ## but only 3 parameters in the (phi, rho) parameterization.
    ## What makes it a "POPAN-lambda" model is that phi and rho are both constant within each group for the
    ## duration of the survey.
    ##
    ## Note also that the code is set up so that per-capita birth and death rates apply within each group. So for
    ## example if the groups are males and females, and if they have the same birth-rate of 0.2, then the
    ## number of males "born" will be 0.2 * (number of current males) each year, and likewise for females.
    ## This observation is only relevant to the POPAN-lambda model which assumes a constant per-capita
    ## population growth rate per year. The more general pent-style models just estimate births per year
    ## as free parameters.
    ##
    ## Additional model blurb below.
    ##
    ## k is the number of capture occasions.
    ##
    ## Nsvec is a vector giving the superpopulation sizes for different groups.  It is a vector, not a list,
    ## and must have the same length as the number of groups, with Nsvec[1] being the Ns for group 1,
    ## Nsvec[2] being the Ns for group 2, etc. The length of Nsvec determines the number of groups (ngp).
    ## Nsvec is treated differently from the other parameters because it is used to determine ngp.
    ##
    ## Each of phi, rho, and p can be varied by group. They all come in lists assumed to be ordered
    ## by group: e.g. phiList[[1]] applies to group 1, phiList[[2]] applies to group 2, etc. Any names given to
    ## these list components are ignored - i.e. renaming components will not reorder the list.
    ## If there is only one element in the list, this is assumed to apply to all groups.
    ## If they are not supplied as a list at all, e.g. pList=c(0.4, 0.5, 0.6), then this will be recast as a list with
    ## first element c(0.4, 0.5, 0.6).

    ## Determine number of groups from the Ns vector:
    ngp <- length(Nsvec)

    ## Ensure phiList, rhoList, and pList have the right number of entries, which should be either 1 or
    ## ngp in each case:
    if(ngp>1){
        for(parname in c("phi", "rho", "p")){
            parListName <- paste0(parname, "List")
            parList <- eval(parse(text=parListName))
            ## If parList is not already in list format, recast it as the first element of a list:
            ## e.g. if pList=c(0.4, 0.5, 0.6) then this will become pList[[1]] = c(0.4, 0.5, 0.6).
            if(!is.list(parList)) parList <- list(parList)
            ## If there is only one element in parList, recycle it for all groups:
            if(length(parList)==1){
                parListNew <- vector("list", ngp)
                for(gp in 1:ngp) parListNew[[gp]] <- parList[[1]]
                assign(x=parListName, value=parListNew)
            }
            else if(length(parList)<ngp) stop(paste(parListName, "should be a list of length 1 or ngp."))
        }
    }

    ## Check the lengths of each list element:
    if(any(lapply(phiList, length)!=1)) stop("phiList should have all elements of length 1")
    if(any(lapply(rhoList, length)!=1)) stop("rhoList should have all elements of length 1")
    if(any(lapply(pList, length)!=k)) stop("pList should have all elements of length k")
    ## Each of phiList, rhoList, and pList are now in the required format of a list with ngp elements.

    ## ----------------------------------------------------------------------------------------------------------------------
    ## Create capture histories for a single group: this uses phiList[[gp]], rhoList[[gp]], and pList[[gp]]:
    caphistGroup.func <- function(gp){
        ## Find pentvec for this group: pentLambda.func is the POPAN-lambda model parametrized by rho.
        pentvec <- pentLambda.func(rho=rhoList[[gp]], phi=phiList[[gp]], k=k)

        ## Find numbers entered on each occasion:
        nentered <- cumsum(rmultinom(1, Nsvec[gp], pentvec))

        ## Create matrices and enter first animals
        alive <- caphists <- matrix(0, nrow=Nsvec[gp], ncol=k)
        alive[1:nentered[1], 1] <- 1

        ## Enter first captures
        ## NB I think there was an error in Robin's code previously: it should be pvec[1], not just pvec below
        ## (recast as pList[[gp]][1] below).
        caphists[1:nentered[1], 1] <- rbinom(nentered[1], 1, pList[[gp]][1])

        ## Loop over remaining occasions
        for (t in 2:k) {
            ## Enter survivors from previous occasion
            alive[which(alive[1:nentered[t - 1], t - 1] == 1), t] <-
                rbinom(sum(alive[1:nentered[t - 1], t - 1]), 1, phiList[[gp]])

            ## Enter new animals
            if(nentered[t]>nentered[t-1]) alive[(nentered[t - 1] + 1):nentered[t], t] <- 1

            ## Find captures
            caphists[which(alive[1:nentered[t], t] == 1), t] <-
                rbinom(sum(alive[1:nentered[t], t]), 1, pList[[gp]][t])
        }

        ## Remove unobserved animals
        caphists <- caphists[rowSums(caphists) > 0, ]
        return(caphists)
    }
    ## ---------------------------------------------------------------------------------------------
    ## Return list of group capture histories:
    lapply(1:ngp, caphistGroup.func)
}


##############################################################

popanLambda.fit.func <- function(dat, k=ncol(dat[[1]]), model=list(
                                                            gp1=c("Ns.1", "rho.1", "phi.1", paste0("p", 1:k, ".1")),
                                                            gp2=c("Ns.2", "rho.1", "phi.1", paste0("p", 1:k, ".2"))),
                                 ## startvec is a named vector: it needs all the entries specified in model, but in any order.
                                 ## If startvec includes parameters that aren't needed for the specified model, they are ignored.
                                 startvec=c(Ns.1=950, Ns.2=1150, rho.1=0.15, phi.1=0.93, rho.2=0.15, phi.2=0.93,
                                            structure(rep(0.45, k), .Names=paste0("p", 1:k, ".1")),  ## start p's for group 1
                                            structure(rep(0.45, k), .Names=paste0("p", 1:k, ".2"))),  ## start p's for group 2
                                 printit=FALSE){

    ## popanLambda.fit.func 11/5/2020
    ##
    ## This is the POPAN-lambda model, but it is parametrized by rho = lambda-phi, which is the per-capita
    ## birth rate (true births plus immigration) for each group.  The reason for switching from lambda to rho
    ## is that rho makes more sense once we have multiple groups: see the blurb for
    ## simGroups.popanLambda.func for more details.
    ##
    ## This function calculates the combined binomial likelihood for the multi-group POPAN-lambda model.
    ## Based on Robin's function twogroupnll, recast by Rachel for multiple groups and generic
    ## parameter structures.
    ##
    ## The data, 'dat', is a list of capture histories by group: e.g. dat[[1]] is a n1 x k matrix of capture histories
    ## for n1 animals belonging to group 1; dat[[2]] is a n2 x k matrix for group 2; etc.
    ## If dat is just a single data-frame, it will be coerced into a list with one element (i.e. just one group).
    ##
    ## Any combination of parameters can be shared between groups by defining the corresponding names
    ## in the model argument. Each group is given parameters in the order Ns, rho, phi, pvec.
    ## Naming is done by parameter.gp, e.g. Ns.1 means the Ns parameter for group 1.
    ## Any parameters that are shared between groups can be specified as such via the model strings:
    ## e.g. if model$gp2=c("Ns.2", "rho.1", ...) this means that Ns.2 is estimated separately from Ns.1, but
    ## rho.2=rho.1 so the two groups share the same rho parameter (within-group per-capita birth-rate).
    ##
    ## The p-vectors are named as poccasion.gp, e.g. p5.2 means capture probability for occasion 5, group 2.
    ## To use a constant capture probability by occasion it would be something like,
    ## model$gp1=c("Ns.1", "rho.1", "phi.1", rep("p1.1", k))

    ## Coerce dat into a list if it isn't already one:
    if(!is.list(dat)) dat <- list(dat)

    ## Check that all elements of dat have the same number of columns, k:
    kcheck <- lapply(dat, ncol)
    if(any(kcheck!=k)) stop("Some elements of the data-list 'dat' have different k or don't match the k supplied.")

    ## Find the number of groups: this is the number of elements in the list 'dat'.
    ngp <- length(dat)

    ## -----------------------------------------------------------------------------------------------------------
    ## Parameter organisation:
    ## -----------------------------------------------------------------------------------------------------------
    ## Create the maximal vector of parameter-names. This is a separate Ns, rho, phi, pvec-t for each group.
    parnamesTemplate <- c("Ns", "rho", "phi", paste0("p", 1:k))
    allparnames <- unlist(lapply(1:ngp, function(gp) paste(parnamesTemplate, gp, sep=".")))
    ## The actual parameters to be estimated are those in unlist(model[1:ngp]).
    ## Create a vector called allparvec that lists all the available parameter-slots in allparnames, and
    ## specifies the name of the parameter to fill that slot. For example, if ngp=2 and k=3,
    ## we might have allparvec as follows: this specifies rho and phi in common between the two groups
    ## and the other parameters separate.
    ## names = Ns.1 rho.1 phi.1 p1.1 p2.1 p3.1 Ns.2 rho.2 phi.2 p1.2 p2.2 p3.2
    ## entries = Ns.1 rho.1 phi.1 p1.1 p2.1 p3.1 Ns.2 rho.1 phi.1 p1.2 p2.2 p3.2
    allparvec <- structure(unlist(model[1:ngp]), .Names=allparnames)

    ## Finally, parsEstvec lists only the parameters that will be estimated, for sending through nlminb:
    ## e.g. in the example above, pars.estvec =c( Ns.1 rho.1 phi.1 p1.1 p2.1 p3.1 Ns.2 p1.2 p2.2 p3.2 )
    ## (all as character strings).
    parsEstvec <- unique(allparvec)
    ## parInds is a vector of length equal to allparnames that specifies which of the estimated parameters
    ## in parsEstvec is used for each parameter in the master-list.  For example if rho.2=rho.1, then
    ## the entry of parInds corresponding to rho.2 would be the same as that corresponding to rho.1
    ## (i.e. index 2). We will populate the whole vector of parameters inside the likelihood by using
    ## allpars = pars[parInds].
    parInds <- match(allparvec, parsEstvec)

    ## Check that all parameters for estimation have been given start values in startvec:
    if(!all(parsEstvec %in% names(startvec))) stop("some parameter start values are missing from startvec.")
    startvec <- startvec[parsEstvec]

    ## -----------------------------------------------------------------------------------------------------
    ## Summarise the data using Robin's summary functions:
    popsumList <- lapply(dat, popansummary.func)
    nhistList <- lapply(popsumList, function(x) attributes(x)$nhist)

    ## -----------------------------------------------------------------------------------------------------
    ## Vectors of lower and upper parameter limits for the nlminb optimization.
    ## -----------------------------------------------------------------------------------------------------
    ## All parameters except for Ns's have lower limit of 0 (set to 1e-6 to keep out of trouble at the boundaries):
    npar <- length(parsEstvec)
    lowervec <- rep(1e-6, npar)
    ## For the Ns parameters, insert lower limit which is the maximum of nhist for all the groups that
    ## possess that parameter.  For example if Ns.1=Ns.2 in the model, then the constraint for Ns.1
    ## is that it must be greater than both nhistList[[1]] and nhistList[[2]].
    for(gp in 1:ngp){
        Nsgp <- paste0("Ns.", gp)
        groupsUsingNsgp <- unique(as.numeric(gsub(pattern="Ns.", replacement="",
                                                  x=allparnames[allparvec==Nsgp])))
        if(length(groupsUsingNsgp)>0)
            lowervec[parsEstvec==Nsgp] <- max(unlist(nhistList[groupsUsingNsgp]))
    }
    ## uppervec is 1 for phi and p parameters, Inf for Ns and rho parameters:
    uppervec <- rep(1-1e-6, npar)
    uppervec[grep("Ns", parsEstvec)] <- Inf
    uppervec[grep("rho", parsEstvec)] <- Inf


    ## -----------------------------------------------------------------------------------------------------
    ## Negative log likelihood function:
    ## -----------------------------------------------------------------------------------------------------
    negloglik.func <- function(pars) {
        ## pars contains only the parameters for estimation.
        ## Populate the vector of all parameters by mapping pars to allpars as follows:
        allpars <- pars[parInds]
        names(allpars) <- allparnames

        ## Inner function to find the neg-log-likelihood contribution for a single group.
        onegroup.func <- function(gp){
            ## Unpack the parameters for this group:
            Ns <- allpars[paste0("Ns.", gp)]
            rho <- allpars[paste0("rho.", gp)]
            phi <- allpars[paste0("phi.", gp)]
            pvec <- allpars[paste0("p", 1:k, ".", gp)]

            ## Find entry proportions from the rho-parametrization of the POPAN-lambda model:
            pentvec <- pentLambda.func(rho, phi, k)

            ## Find psi, chi, and ptheta for this group:
            psivec <- psi.func(pentvec, rep(phi, k-1), pvec, k)  ## phi is constant in the POPAN-lambda model
            chivec <- chi.func(rep(phi, k-1), pvec, k)
            ptheta <- ptheta.func(pentvec, pvec, chivec)

            ## Number of capture histories in this group:
            popsum <- popsumList[[gp]]
            nhist <- nhistList[[gp]]

            ## Return the negative log-likelihood contribution for this group:
            -sum(
                 ## Binomial coefficients
                 ## Could use lchoose(Ns, nhist) instead of the next three lines; it gives slightly different answers
                 ## (<=1% difference in parameter estimates over the course of the optimization) even though
                 ## it should be identical.
                 lgamma(Ns + 1),
                 -lgamma(nhist + 1),
                 -lgamma(Ns - nhist + 1),
                 ## Group capture histories including undetected
                 (Ns - nhist) * log(1 - ptheta),
                 popsum$first.tab * log(psivec),
                 popsum$caps * log(pvec),
                 popsum$non.caps * log(1 - pvec),
                 popsum$survives[-k] * log(phi),
                 popsum$last.tab * log(chivec))
        }
        ## Overall negative-log-likelihood is the sum across groups:
        nll <- sum(sapply(1:ngp, onegroup.func))
        if(printit){
            for (i in 1:length(pars)){
                cat(names(pars)[i], ": ", round(pars[i], 3), ", ", sep = "")
            }
            cat("NLL:", nll, "\n")
        }
        return(nll)
    }

    ## --------------------------------------------------------------------------------------------------
    ## Fit popan model and return:
    ## -----------------------------------------------------------------------------------------------------
    nlminb(
        start = startvec,
        objective = negloglik.func,
        lower = lowervec,
        upper = uppervec)


}

##############################################################

popanLambda.wrap <- function(Nsim=100, k=11,  # k=number of capture occasions
                             Nsvec=c(1000, 1200), phiList=0.85, rhoList=0.18, pList=list(rep(0.3, k), rep(0.2, k)), #true vals
                             model=list( ## Model to fit: common phi and rho, different Ns and p by group; p const over time
                                 gp1=c("Ns.1", "rho.1", "phi.1", rep("p1.1", k)),
                                 gp2=c("Ns.2", "rho.1", "phi.1", rep("p1.2", k))),
                             ## start from the generating values: the optimization needs good start-points.
                             startvec=c(Ns.1=1000, Ns.2=1200, rho.1=0.18, phi.1=0.85, p1.1=0.3, p1.2=0.2)){
    ## popanLambda.wrap 12/5/2020
    ## Repeatedly simulate from simGroups.popanLambda.func and fit the model.
    ## Uses the (phi,rho) parametrization of POPAN-lambda: see the blurb for simGroups.popanLambda.func.
    ##
    ## This function doesn't do any of the usual error-catching, as Ben might want to do his own thing.

    onesim.func <- function(sim){
        ## Single simulation:
        if(sim%%10==0) cat("Simulation", sim, "\n")
        simdat <- simGroups.popanLambda.func(k=k,
                                             Nsvec=Nsvec, phiList=phiList, rhoList=rhoList, pList=pList)

        ## Fit popan-Lambda model specified:
        popanLambda.fit.func(simdat, k=k, model=model, startvec=startvec)$par

    }
    ## Run Nsim simulations and return the data-frame of estimates:
    as.data.frame(t(sapply(1:Nsim, onesim.func)))

}

##############################################################

immigrationElNino.func <- function(bpars=c(b1=0.18, b2=0.1)){
    ## immigrationElNino.func 13/5/2020
    ## Example function to show how we can incorporate covariates into pent parameters.
    ## This function takes a vector of parameters, bpars for 'birth parameters', and returns a vector
    ## rhovec of birth rates per year.
    ## The specific example given here is that there is a constant background birth-rate per year (b1),
    ## consistent with the POPAN-lambda model, but in ElNino years there is a boost of b2 due to immigration.
    ## The ElNino years are input here as 2014-2015, which will generate an influx to the population for
    ## 2015-2016: this is just Rachel's guess based on Edy's emails.
    ##
    ## So rho=b1 for all years 2009, 2010, 2011, 2012, 2013, 2016, 2017, 2018,
    ## and rho=b1+b2 for 2014, 2015.
    ## Note that rhovec only has length k-1: the final entry for 2018 controls new animals that enter in 2019
    ## as a per-capita proportion of the 2018 animals. Thus for a population size of N_2018 in 2018, the
    ## expected number of animals entering the population for the first time in 2019 is rho_2018*N_2018.
    ##
    ## If the population is being modelled in groups (e.g. males and females), this per-capita proportion works
    ## within each group: so the number of new males is proportional to the current number of males, and
    ## the number of new females is proportional to the current number of females. Although this sounds odd,
    ## it does fit a model with a single birth-rate dependent only on the number of mothers.  For example,
    ## suppose the sex ratio is m = P(male), so expected population sizes satisfy Nm/(Nm+Nf)=m,
    ## giving m*Nf = (1-m)*Nm.
    ## If new births = alpha*Nf, then new males = alpha*m*Nf =alpha*(1-m)*Nm, so per-capita birth-rate in
    ## the males group is alpha*(1-m).
    ## Likewise, new females=alpha*(1-m)*Nf, so per-capita birth-rate in the females group is also alpha*(1-m).
    ##
    ## NOTE: there is nothing in this setup that allows immigrants to be treated differently once they've
    ## immigrated.  If the immigration is only temporary and these same individuals leave the population
    ## after 2016, then a two-point mixture might be needed for the subsequent survival probabilities, phi.
    ## The correct mixture proportions should be calculable from b1 and b2 so they wouldn't need to be
    ## estimated.

    ## Create the ElNino data-frame (noting this is just Rachel's guesswork):
    ## usually the covariate data would be externally-defined and supplied to this function through arguments.
    elnino.dat <- data.frame(year=2009:2018, elnino=rep(0, 10))
    elnino.dat$elnino[elnino.dat$year %in% c(2015, 2016)] <- 1

    ## The required rho-vector is b1 + elnino * b2:
    rhovec <- bpars[1] + bpars[2]*elnino.dat$elnino
    return(rhovec)
}

###############################################################

simGroups.popanGeneral.func <- function(k=11,
                                        Nsvec=c(1000, 1200),
                                        phiList=list(rep(0.85, k-1), rep(0.85, k-1)),  ## Can enter as phiList=list(0.85) if const.
                                        ## Alternatively, install a function "survfunc" for phi, like "birthfunc" replaces rho below.
                                        birthfunc=immigrationElNino.func, ## Custom covariate function for births/immigrants
                                        birthparList=list(bpar1=c(0.18, 0.1), bpar2=c(0.18, 0.1)),
                                        ## birthparList specifies inputs to birthfunc for each group
                                        pList= list(rep(0.3, k), rep(0.2, k))
                                        ## If p is just one number for a group, it will be taken as constant for the group.
                                        ){
    ## simGroups.popanGeneral.func 13/5/2020
    ## Simulate from the POPAN-general model for multiple groups (e.g. males and females), in which
    ## the birth-curve is input through a function birthfunc that can take covariates.
    ## Currently only birthfunc (replacing rhoList) has been made into a custom-covariate function,
    ## but phiList and pList could also be replaced by functions if desired.
    ##
    ## This function is based on template function simGroups.popanLambda.func. To use this function to
    ## replicate the POPAN-lambda model, use birthfunc=function(b1) rep(b1, k-1) and
    ## (say) birthparList=list(0.18, 0.18) : i.e. constant per-capita birth-rate within each group.
    ##
    ## Note that the code is set up so that per-capita birth and death rates apply within each group. So for
    ## example if the groups are males and females, and if they have the same birth-rate of 0.2, then the
    ## number of males "born" will be 0.2 * (number of current males) each year, and likewise for females.
    ## See the blurb for immigrationElNino.func to see why this makes sense even if the number of new males
    ## is completely reliant on the number of females at the previous step, rather than the number of males,
    ## if we have a population with non-equal sex-ratio.
    ##
    ## Additional model blurb below.
    ##
    ## k is the number of capture occasions.
    ##
    ## Nsvec is a vector giving the superpopulation sizes for different groups.  It is a vector, not a list,
    ## and must have the same length as the number of groups, with Nsvec[1] being the Ns for group 1,
    ## Nsvec[2] being the Ns for group 2, etc. The length of Nsvec determines the number of groups (ngp).
    ## Nsvec is treated differently from the other parameters because it is used to determine ngp.
    ##
    ## Each of phi, birthpar, and p can be varied by group. They all come in lists assumed to be ordered
    ## by group: e.g. phiList[[1]] applies to group 1, phiList[[2]] applies to group 2, etc. Any names given to
    ## these list components are ignored - i.e. renaming components will not reorder the list.
    ## If there is only one element in the list, this is assumed to apply to all groups.
    ## If they are not supplied as a list at all, e.g. pList=c(0.4, 0.5, 0.6), then this will be recast as a list with
    ## first element c(0.4, 0.5, 0.6).

    ## Determine number of groups from the Ns vector:
    ngp <- length(Nsvec)

    ## Ensure phiList, birthparList, and pList have the right number of entries, which should be either 1 or
    ## ngp in each case:
    if(ngp>1){
        for(parname in c("phi", "birthpar", "p")){
            parListName <- paste0(parname, "List")
            parList <- eval(parse(text=parListName))
            ## If parList is not already in list format, recast it as the first element of a list:
            ## e.g. if pList=c(0.4, 0.5, 0.6) then this will become pList[[1]] = c(0.4, 0.5, 0.6).
            if(!is.list(parList)) parList <- list(parList)
            ## If there is only one element in parList, recycle it for all groups:
            if(length(parList)==1){
                parListNew <- vector("list", ngp)
                for(gp in 1:ngp) parListNew[[gp]] <- parList[[1]]
                assign(x=parListName, value=parListNew)
            }
            else if(length(parList)<ngp) stop(paste(parListName, "should be a list of length 1 or ngp."))
        }
    }

    ## Within phiList, all elements should have either 1 or k-1 entries: if 1 entry we take phi as constant
    ## and recycle it to length k-1.
    phiList <- lapply(phiList, function(x) if(length(x)==1) return(rep(x, k-1)) else return(x))
    ## Within pList, all elements should have either 1 or k entries: if 1 entry take p as constant and recycle it.
    pList <- lapply(pList, function(x) if(length(x)==1) return(rep(x, k)) else return(x))

    ## ----------------------------------------------------------------------------------------------------------------------
    ## Convert birthparList to rhoList using the function birthfunc supplied:
    rhoList <- lapply(birthparList, birthfunc)

    ## -------------------------------------------------------------------------------------------------------------
    ## Check the lengths of each list element: all phi and rho vectors should have length k-1,
    ## p vectors should have length k:
    if(any(lapply(phiList, length)!=k-1)) stop("phiList should have all elements of length k-1")
    if(any(lapply(rhoList, length)!=k-1)) stop("rhoList should have all elements of length k-1")
    if(any(lapply(pList, length)!=k)) stop("pList should have all elements of length k")
    ## Each of phiList, rhoList, and pList are now in the required format of a list with ngp elements.

    ## ----------------------------------------------------------------------------------------------------------------------
    ## Create capture histories for a single group: this uses phiList[[gp]], rhoList[[gp]], and pList[[gp]]:
    caphistGroup.func <- function(gp){
        ## Find pentvec for this group: pentGeneral.func takes a general phivec and rhovec and returns pentvec:
        pentvec <- pentGeneral.func(rhovec=rhoList[[gp]], phivec=phiList[[gp]], k=k)

        ## Find numbers entered on each occasion:
        nentered <- cumsum(rmultinom(1, Nsvec[gp], pentvec))

        ## Create matrices and enter first animals
        alive <- caphists <- matrix(0, nrow=Nsvec[gp], ncol=k)
        alive[1:nentered[1], 1] <- 1

        ## Enter first captures
        caphists[1:nentered[1], 1] <- rbinom(nentered[1], 1, pList[[gp]][1])

        ## Loop over remaining occasions
        for (t in 2:k) {
            ## Enter survivors from previous occasion: use phivec[t-1] for survival from time t-1 to time t:
            alive[which(alive[1:nentered[t - 1], t - 1] == 1), t] <-
                rbinom(sum(alive[1:nentered[t - 1], t - 1]), 1, phiList[[gp]][t-1])

            ## Enter new animals
            if(nentered[t]>nentered[t-1]) alive[(nentered[t - 1] + 1):nentered[t], t] <- 1

            ## Find captures
            caphists[which(alive[1:nentered[t], t] == 1), t] <-
                rbinom(sum(alive[1:nentered[t], t]), 1, pList[[gp]][t])
        }

        ## Remove unobserved animals
        caphists <- caphists[rowSums(caphists) > 0, ]
        return(caphists)
    }
    ## ---------------------------------------------------------------------------------------------
    ## Return list of group capture histories:
    lapply(1:ngp, caphistGroup.func)
}


##############################################################

popanGeneral.fit.func <- function(dat, k=ncol(dat[[1]]), birthfunc = immigrationElNino.func,
                                  model=list(
                                      gp1=c("Ns.1", "b1.1", "b2.1", rep("phi1.1", k-1), paste0("p", 1:k, ".1")),
                                      gp2=c("Ns.2", "b1.1", "b2.1", rep("phi1.1", k-1), paste0("p", 1:k, ".2"))),
                                  ## Note: parameters must be entered into model in the order Ns, b, phi, p within each group.
                                  ## The default model above allows different group sizes, time-varying capture probabilities
                                  ## that differ by group, a single phi for all times and both groups, and a common per-capita
                                  ## birth rate within each group. The parameter naming convention is occasion.group,
                                  ## e.g. phi5.1 is the phi for occasion 5 in group 1.
                                  ##
                                  ## startvec is a named vector: it needs all the entries specified in model, but in any order.
                                  ## If startvec includes parameters that aren't needed for the specified model, they are ignored.
                                  startvec=c(Ns.1=950, Ns.2=1150, b1.1=0.15, b2.1=0.1, phi1.1=0.93,
                                             structure(rep(0.45, k), .Names=paste0("p", 1:k, ".1")),  ## start p's for group 1
                                             structure(rep(0.45, k), .Names=paste0("p", 1:k, ".2"))),  ## start p's for group 2
                                  printit=FALSE){

    ## popanGeneral.fit.func 13/5/2020
    ## Fit the POPAN-general model for multiple groups (e.g. males and females).
    ## The birth-curve is input through a function birthfunc that can take covariates. There can be any number
    ## of parameters in birthfunc, but they should be labelled b1, b2, b3, etc.  Any parameter labelled b* will be
    ## assumed to be an argument for birthfunc. The parameters are fed to birthfunc in a single vector, i.e.
    ## birthfunc has argument bpars which is a vector of parameter inputs, e.g. bpars=c(b1, b2).
    ##
    ## Currently only birthfunc (replacing rho) has been made into a custom-covariate function,
    ## but phi and p could also be replaced by functions if desired.
    ##
    ## To use this function to replicate the POPAN-lambda model, use birthfunc=function(b1) rep(b1, k-1).
    ##
    ## See the blurb for simGroups.popanGeneral.func for more details on the model and code structure.
    ##
    ## The data, 'dat', is a list of capture histories by group: e.g. dat[[1]] is a n1 x k matrix of capture histories
    ## for n1 animals belonging to group 1; dat[[2]] is a n2 x k matrix for group 2; etc.
    ## If dat is just a single data-frame, it will be coerced into a list with one element (i.e. just one group).
    ##
    ## Any combination of parameters can be shared between groups by defining the corresponding names
    ## in the model argument. Each group is given parameters in the order Ns, b, phi, pvec.
    ## Naming is done by parameter.gp, e.g. Ns.1 means the Ns parameter for group 1.
    ## Any parameters that are shared between groups can be specified as such via the model strings:
    ## e.g. if model$gp2=c("Ns.2", "b1.1", ...) this means that Ns.2 is estimated separately from Ns.1, but
    ## b1.2=b1.1 so the two groups share the same b1 parameter (the first argument of function birthfunc).
    ##
    ## Vectors are named as occasion.gp, e.g. p5.2 means capture probability for occasion 5, group 2.
    ## For constant phi and constant capture probability by occasion, input something like:
    ## model$gp1=c("Ns.1", "b1.1", "b2.1", rep("phi1.1", k-1), rep("p1.1", k))

    ## Coerce dat into a list if it isn't already one:
    if(!is.list(dat)) dat <- list(dat)

    ## Check that all elements of dat have the same number of columns, k:
    kcheck <- lapply(dat, ncol)
    if(any(kcheck!=k)) stop("Some elements of the data-list 'dat' have different k or don't match the k supplied.")

    ## Find the number of groups: this is the number of elements in the list 'dat'.
    ngp <- length(dat)

    ## -----------------------------------------------------------------------------------------------------------
    ## Parameter organisation:
    ## -----------------------------------------------------------------------------------------------------------
    ## Create the maximal vector of parameter-names. This is a separate Ns, bvec, phivec, pvec for each group.
    ## Detect how many b-parameters are required for birthfunc using the formals() function:
    ## this assumes that birthfunc takes a single vector argument with a default of the correct length
    ## specified, and we are aiming to find how long that default is:
    nbpar <- length(eval(formals(birthfunc)[[1]]))
    parnamesTemplate <- c("Ns", paste0("b", 1:nbpar), paste0("phi", 1:(k-1)), paste0("p", 1:k))
    allparnames <- unlist(lapply(1:ngp, function(gp) paste(parnamesTemplate, gp, sep=".")))

    ## The actual parameters to be estimated are those in unlist(model[1:ngp]).
    ## Create a vector called allparvec that lists all the available parameter-slots in allparnames, and
    ## specifies the name of the parameter to fill that slot. For example, if ngp=2 and k=3,
    ## we might have allparvec as follows: this specifies b's and phi's in common between the two groups
    ## and the other parameters separate.
    ## names = Ns.1 b1.1 b2.1 phi1.1 phi2.1 p1.1 p2.1 p3.1 Ns.2 b1.2 b2.2 phi1.2 phi2.2 p1.2 p2.2 p3.2
    ## entries = Ns.1 b1.1 b2.1 phi1.1 phi1.1 p1.1 p2.1 p3.1 Ns.2 b1.1 b2.1 phi1.1 phi1.1 p1.2 p2.2 p3.2
    allparvec <- structure(unlist(model[1:ngp]), .Names=allparnames)

    ## Finally, parsEstvec lists only the parameters that will be estimated, for sending through nlminb:
    ## e.g. in the example above, pars.estvec =c( Ns.1 b1.1 b2.1 phi1.1 p1.1 p2.1 p3.1 Ns.2 p1.2 p2.2 p3.2 )
    ## (all as character strings).
    parsEstvec <- unique(allparvec)
    ## parInds is a vector of length equal to allparnames that specifies which of the estimated parameters
    ## in parsEstvec is used for each parameter in the master-list.  For example if b1.2=b1.1, then
    ## the entry of parInds corresponding to b1.2 would be the same as that corresponding to b1.1
    ## (i.e. index 2). We will populate the whole vector of parameters inside the likelihood by using
    ## allpars = pars[parInds].
    parInds <- match(allparvec, parsEstvec)

    ## Check that all parameters for estimation have been given start values in startvec:
    if(!all(parsEstvec %in% names(startvec))) stop("some parameter start values are missing from startvec.")
    startvec <- startvec[parsEstvec]

    ## -----------------------------------------------------------------------------------------------------
    ## Summarise the data using Robin's summary functions:
    popsumList <- lapply(dat, popansummary.func)
    nhistList <- lapply(popsumList, function(x) attributes(x)$nhist)

    ## -----------------------------------------------------------------------------------------------------
    ## Vectors of lower and upper parameter limits for the nlminb optimization.
    ## -----------------------------------------------------------------------------------------------------
    ## Assume the bpars are restricted between 0 and infinity (or edit below if not).
    ## All parameters except for Ns's have lower limit of 0 (set to 1e-6 to keep out of trouble at the boundaries):
    npar <- length(parsEstvec)
    lowervec <- rep(1e-6, npar)
    lowervec[grep("b", parsEstvec)] <- -Inf
    ## For the Ns parameters, insert lower limit which is the maximum of nhist for all the groups that
    ## possess that parameter.  For example if Ns.1=Ns.2 in the model, then the constraint for Ns.1
    ## is that it must be greater than both nhistList[[1]] and nhistList[[2]].
    for(gp in 1:ngp){
        Nsgp <- paste0("Ns.", gp)
        groupsUsingNsgp <- unique(as.numeric(gsub(pattern="Ns.", replacement="",
                                                  x=allparnames[allparvec==Nsgp])))
        if(length(groupsUsingNsgp)>0)
            lowervec[parsEstvec==Nsgp] <- max(unlist(nhistList[groupsUsingNsgp]))
    }
    ## uppervec is 1 for phi and p parameters, Inf for Ns and bpar parameters:
    uppervec <- rep(1-1e-6, npar)
    uppervec[grep("Ns", parsEstvec)] <- Inf
    uppervec[grep("b", parsEstvec)] <- Inf

    ## -----------------------------------------------------------------------------------------------------
    ## Negative log likelihood function:
    ## -----------------------------------------------------------------------------------------------------
    negloglik.func <- function(pars) {
        ## pars contains only the parameters for estimation.
        ## Populate the vector of all parameters by mapping pars to allpars as follows:
        allpars <- pars[parInds]
        names(allpars) <- allparnames

        ## Inner function to find the neg-log-likelihood contribution for a single group.
        onegroup.func <- function(gp){
            ## Unpack the parameters for this group:
            Ns <- allpars[paste0("Ns.", gp)]
            bpars <- allpars[paste0("b", 1:nbpar, ".", gp)]  ## nbpars of these
            phivec <- allpars[paste0("phi", 1:(k-1), ".", gp)]  ## k-1 of these
            pvec <- allpars[paste0("p", 1:k, ".", gp)]  ## k of these

            ## Use the birthfunc supplied to convert bpars into rhovec:
            rhovec <- birthfunc(bpars)

            ## Find entry proportions from the POPAN-general function:
            pentvec <- pentGeneral.func(rhovec, phivec, k)

            ## Find psi, chi, and ptheta for this group:
            psivec <- psi.func(pentvec, phivec, pvec, k)
            chivec <- chi.func(phivec, pvec, k)
            ptheta <- ptheta.func(pentvec, pvec, chivec)

            ## Number of capture histories in this group:
            popsum <- popsumList[[gp]]
            nhist <- nhistList[[gp]]

            ## Return the negative log-likelihood contribution for this group:
            -sum(
                 ## Binomial coefficients
                 ## Could use lchoose(Ns, nhist) instead of the next three lines; it gives slightly different answers
                 ## (<=1% difference in parameter estimates over the course of the optimization) though
                 ## it should be identical.
                 lgamma(Ns + 1),
                 -lgamma(nhist + 1),
                 -lgamma(Ns - nhist + 1),
                 ## Group capture histories including undetected
                 (Ns - nhist) * log(1 - ptheta),
                 popsum$first.tab * log(psivec),
                 popsum$caps * log(pvec),
                 popsum$non.caps * log(1 - pvec),
                 popsum$survives[-k] * log(phivec),
                 popsum$last.tab * log(chivec))
        }
        ## Overall negative-log-likelihood is the sum across groups:
        nll <- sum(sapply(1:ngp, onegroup.func))
        if(printit){
            for (i in 1:length(pars)){
                cat(names(pars)[i], ": ", round(pars[i], 3), ", ", sep = "")
            }
            cat("NLL:", nll, "\n")
        }
        return(nll)
    }


    ## --------------------------------------------------------------------------------------------------
    ## Fit popan model and return:
    ## -----------------------------------------------------------------------------------------------------
    nlminb(
        start = startvec,
        objective = negloglik.func,
        lower = lowervec,
        upper = uppervec,
        control=list(iter.max=1000, eval.max=5000))


}

##############################################################

popanGeneral.wrap <- function(Nsim=100, k=11,  # k=number of capture occasions
                              birthfunc = immigrationElNino.func,
                              Nsvec=c(1000, 1200),
                              phiList=0.85,  ## specifying one phi means it will be cycled for all times and groups
                              birthparList=list(bpar1 = c(0.18, 0.1)),
                              ## just one list element means these parameters will be used for both groups
                              pList= list(rep(0.3, 11), rep(0.2, 11)),  ## two list elements => different p's for the groups
                              ##
                              ## Model to fit: common birth and survival, different Ns and p by group; phi & p const over time
                              model=list(
                                  gp1=c("Ns.1", "b1.1", "b2.1", rep("phi1.1", 10), rep("p1.1", 11)),
                                  gp2=c("Ns.2", "b1.1", "b2.1", rep("phi1.1", 10), rep("p1.2", 11))),
                              ## start from the generating values: the optimization needs good start-points.
                              startvec=c(Ns.1=1000, Ns.2=1200, b1.1=0.18, b2.1=0.1, phi1.1=0.85, p1.1=0.3, p1.2=0.2))
{
    ## popanGeneral.wrap 13/5/2020
    ## Repeatedly simulate from simGroups.popanGeneral.func and fit the model.
    ##
    ## This function doesn't do any of the usual error-catching, as Ben might want to do his own thing.

    onesim.func <- function(sim){
        ## Single simulation:
        if(sim%%10==0) cat("Simulation", sim, "\n")
        simdat <- simGroups.popanGeneral.func(k=k, birthfunc = birthfunc,
                                              Nsvec=Nsvec, phiList=phiList, birthparList=birthparList, pList=pList)

        ## Fit popan-General model specified:
        popanGeneral.fit.func(simdat, k=k, birthfunc = birthfunc, model=model, startvec=startvec)$par

    }
    ## Run Nsim simulations and return the data-frame of estimates:
    as.data.frame(t(sapply(1:Nsim, onesim.func)))

}

##############################################################
