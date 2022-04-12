pentGeneral.transience.func <- function(rhovec, phivec, ptrvec, k){
    ## pentGeneral.transience.func does the same thing as above, but
    ## incorporates transience. On occasion t, the proportion of new
    ## entrants that are transients is given by ptrvec[t].

    gammavec <- numeric(k - 2)
    for (t in 1:(k - 2)){
        gammavec[t] <- phivec[t] + (1 - ptrvec[t + 1])*rhovec[t]
    }
  
    pentvec <- numeric(k)
    pentvec[1] <- 1
    pentvec[2] <- rhovec[1]*(1 - ptrvec[1])
    for (t in 2:(k - 1)){
        pentvec[t + 1] <- rhovec[t]*prod(gammavec[1:(t - 1)])*(1 - ptrvec[1])
    }
    pentvec / sum(pentvec)
}

###############################################################

popanGeneral.covs.fit.func.transience <- function(dat, k=ncol(dat[[1]]), birthfunc = immigrationElNino.func, phifunc, pfunc,
                                       ptrfunc,
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
                                       transience = TRUE,
                                       startvec=c(Ns.1=950, Ns.2=1150, b1.1=0.15, b2.1=0.1, phi1.1=0.93,
                                                  structure(rep(0.45, k), .Names=paste0("p", 1:k, ".1")),  ## start p's for group 1
                                                  structure(rep(0.45, k), .Names=paste0("p", 1:k, ".2"))),
                                       lowervec = NULL, uppervec = NULL,  ## start p's for group 2
                                       printit=FALSE){
    
    ## popanGeneral.fit.func
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
    nphipar <- length(eval(formals(phifunc)[[1]]))
    nppar <- length(eval(formals(pfunc)[[1]]))
    nptrpar <- length(eval(formals(ptrfunc)[[1]]))
    parnamesTemplate <- c("Ns", paste0("b", 1:nbpar), paste0("phi", 1:nphipar),
                          paste0("p", 1:nppar), paste0("ptr", 1:nptrpar)[transience])
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
    startvec.save <- startvec
    startvec <- startvec[parsEstvec]

    ## -----------------------------------------------------------------------------------------------------
    ## Summarise the data using Robin's summary functions:
    popsumList <- lapply(dat, popansummary.func)
    nhistList <- lapply(popsumList, function(x) attributes(x)$nhist)
    first.obsList <- lapply(popsumList, function(x) attributes(x)$first.obs)
    last.obsList <- lapply(popsumList, function(x) attributes(x)$last.obs)
    det.datList <- lapply(popsumList, function(x) attributes(x)$det.dat)
    ## -----------------------------------------------------------------------------------------------------
    ## Vectors of lower and upper parameter limits for the nlminb optimization.
    ## -----------------------------------------------------------------------------------------------------
    ## Assume the bpars are restricted between 0 and infinity (or edit below if not).
    ## All parameters except for Ns's have lower limit of 0 (set to 1e-6 to keep out of trouble at the boundaries):
    npar <- length(parsEstvec)
    if (is.null(lowervec)){
        lowervec <- rep(1e-6, npar)
        lowervec[grep("b", parsEstvec)] <- -Inf
        lowervec[grep("p", parsEstvec)] <- -Inf
        lowervec[grep("phi", parsEstvec)] <- -Inf
        lowervec[grep("ptr", parsEstvec)] <- -Inf
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
    } else {
        lowervec <- lowervec[parsEstvec]
    }
    ## uppervec is 1 for phi and p parameters, Inf for Ns, bpar, and phipar parameters:
    if (is.null(uppervec)){
        uppervec <- rep(1-1e-6, npar)
        uppervec[grep("Ns", parsEstvec)] <- Inf
        uppervec[grep("b", parsEstvec)] <- Inf
        uppervec[grep("p", parsEstvec)] <- Inf
        uppervec[grep("phi", parsEstvec)] <- Inf
        uppervec[grep("ptr", parsEstvec)] <- Inf
    } else {
        uppervec <- uppervec[parsEstvec]
    }
    ## -----------------------------------------------------------------------------------------------------
    ## Negative log likelihood function:
    ## -----------------------------------------------------------------------------------------------------
    negloglik.func <- function(pars, out = "nll") {

        ## pars contains only the parameters for estimation.
        ## Populate the vector of all parameters by mapping pars to allpars as follows:
        allpars <- pars[parInds]
        names(allpars) <- allparnames

        ## Inner function to find the neg-log-likelihood contribution for a single group.
        onegroup.func <- function(gp, out){
            ## Unpack the parameters for this group:
            Ns <- allpars[paste0("Ns.", gp)]
            bpars <- allpars[paste0("b", 1:nbpar, ".", gp)]  ## nbpars of these
            phipars <- allpars[paste0("phi", 1:nphipar, ".", gp)]  ## nphipar of these
            ppars <- allpars[paste0("p", 1:nppar, ".", gp)]  ## nppar of these
            if (transience){
                ptrpars <- allpars[paste0("ptr", 1:nptrpar, ".", gp)]
            }
            ## Use the birthfunc supplied to convert bpars into rhovec:
            rhovec <- birthfunc(bpars)
            ## Use phifunc supplied to convert phipars into phivec.
            phivec <- phifunc(phipars)
            ## Use phifunc supplied to convert phipars into phivec.
            pvec <- pfunc(ppars)
            ## Use ptrfunc supplied to convert ptrpars to ptrvec.
            if (transience){
                ptrvec <- ptrfunc(ptrpars)
            } else {
                ptrvec <- rep(0, k - 1)
            }
            ## Final element of ptrvec doesn't affect the likelihood.
            ptrvec.full <- c(ptrvec, 0)
            ## Find entry proportions from the POPAN-general function:
            pentvec <- pentGeneral.transience.func(rhovec, phivec, ptrvec, k)
            ## Number of capture histories in this group:
            popsum <- popsumList[[gp]]
            nhist <- nhistList[[gp]]
            first.obs <- first.obsList[[gp]]
            last.obs <- last.obsList[[gp]]
            det.dat <- det.datList[[gp]]
            if (Ns < nhist) return(NA)
            
            liktype <- "cpp"
            if (liktype == "rachel"){
                ## Find the chi parameters.  chivec[t] = P(never seen after occasion t | alive at t).
                pentvec <- pentGeneral.func(rhovec, phivec, k)
                ## Simultaneously with chivec, create a vector of length k-1 called psurvive.gap, such that
                ## psurvive.gap[i] = P(survive from survey i to survey i+1 | alive at survey i)
                ## and it is the product of phivec probabilities spanning the calendar years from survey i to survey i+1.
                chivec <- numeric(k)
                chivec[k] <- 1
                psurvive.gap <- numeric(k-1)
                phi.min.index <- length(phivec) + 1
                for(i in (k-1):1) {
                        chivec[i] <- 1 - phivec[i] + phivec[i] * (1-pvec[i+1]) * chivec[i+1]
                }
                log.phivec <- log(phivec)

                ## First create the Binomial portion of the likelihood, and account for all animals never seen.
                ## p.unseen = pent1 (1-p1) chi_1 + pent2 (1-p2) chi_2 + ... + pent.k (1-pk) chi_k
                log.p.unseen <- log(sum(pentvec * (1-pvec) * chivec))
                ## The N-nhist animals never seen give likelihood contribution proportional to
                ## N! / (N-nhist)! p.unseen^(N-nhist)
                ## so the analytic continuation of the negative log likelihood contribution is
                ## -lgamma(N+1) + lgamma(N-nhist+1) - (N-nhist) log(p.unseen):
                nllike <- -lgamma(Ns + 1) + lgamma(Ns - nhist + 1) - (Ns - nhist) * log.p.unseen
                ## Constant so that we get the same log-likelihood as Robin's code.
                nllike <- nllike + lfactorial(nhist)
                ## Now for the animals that were seen.
                ## Go through the data frame according to the occasion of first sighting,
                ## finding capture history probabilities:
                for(f in 1:k) if(any(first.obs==f)){
                        ## (f for first-obs)
                        ## Extract the data records with first sighting on occasion f:
                        dat.f <- det.dat[first.obs==f, ,drop=F]
                        last.dat.f <- last.obs[first.obs==f]
                        n.f <- length(last.dat.f)

                        ## Calculate the psi parameters for this value of first-obs.
                        ## For a particular f, psif.vec is a vector of length f, where
                        ## psif.vec[i] = P(not seen on occasions i, ..., f-1 and survives to occasion f | alive at i)
                        psif.vec <- numeric(f)
                        psif.vec[f] <- 1
                        if(f > 1) for( i in (f-1):1) psif.vec[i] = (1 - pvec[i]) * phivec[i] * psif.vec[i+1]

                        ## All histories with first observation at time f have the same probability up to time f-1:
                        ## prob.to.f.m.1 = (pent1 * psif.1 + ... + pent.f * psif.f)
                        prob.to.f.m.1 <- sum(pentvec[1:f] * psif.vec)
                        log.prob.to.f.m.1 <- log(prob.to.f.m.1)

                        ## Find the capture history probabilities.  For the capture history
                        ##               F       L
                        ## (0  ...  0  1  ...  1  0 ... 0)
                        ##
                        ## where the record is delta_i for i=F, ..., L, the probability is:
                        ##
                        ## (i) if L=F,    prob.to.f.m.1 * p_L * chi_L
                        ## (ii) if L >= F+1,
                        ##         prob.to.f.m.1 * p_L * chi_L *
                        ##                prod_{i=F}^{L-1} [ p_i^delta_i  (1-p_i)^{1-delta_i} * phivec_i ]

                        ## First: ALL animals contribute prob.to.f.m.1 * pvec[L] * chivec[L], so add this in for all animals:
                        #nllike <- nllike - n.f * (log.prob.to.f.m.1)  - sum(log(pvec[last.dat.f])) - sum(log(chivec[last.dat.f]))

                        ## Now add in the extra probabilities of records from occasions f to L-1, for animals
                        ## that have L >= F+1:
                        for(hst in 1:n.f) {
                                last.hst <- last.dat.f[hst]
                                ## datvec.hst is the string of 0s and 1s for the single history hst, from the occasion f of its
                                ## first sighting to the occasion L-1 immediately previous to its last sighting.
                                datvec.hst <- dat.f[hst, f:last.hst]
                                p.hst <- pvec[f:last.hst]
                                log.phivec.hst <- log.phivec[f:(last.hst - 1)]
                                ## The contribution to the log-likelihood from datvec.hst is
                                ## log { prod_{i=F}^{L-1} [ p_i^delta_i  (1-p_i)^{1-delta_i} * phivec_i ] }
                                ## = sum_{i=F}^{L-1} {  delta_i*log(p_i) + (1-delta_i)*log(1-p_i) + log(phivec_i)}

                                ## first bit
                                nllike <- nllike - log.prob.to.f.m.1
                                ## main bit
                                nllike <- nllike - sum(datvec.hst * log(p.hst) + (1 - datvec.hst) * log(1 - p.hst))
                                ## Survival between first and last detections if we have more than one.
                                if (sum(datvec.hst) > 1){
                                    nllike <- nllike - sum(log.phivec.hst)
                                }
                                ## after bit
                                nllike  <- nllike - log(chivec[last.hst])
                                        #nllike <- nllike - sum(loglik.hst)
                        }

                }  ## End of f (likelihood for all animals with first sighting at occasion f).

                #if(is.na(nllike)) nllike <- 20000*(max(abs(pars), na.rm=T)+1)
                #if(is.infinite(nllike)) nllike <- 30000*(max(abs(pars), na.rm=T)+1)

                ## Print current values if required:
                ## print(data.frame(c(allparvec, "nllike"), c(allvalues, nllike)));cat("\n")
            }  else if (liktype == "ben"){
            ## Probability of not being seen after occasion t, given alive at occasion t.
            chivec.trans <- rep(1, k)
            chivec.resid <- numeric(k)
            chivec.resid[k] <- 1
            for (i in (k - 1):1){
                chivec.resid[i] <- 1 - phivec[i] + phivec[i]*(1 - pvec[i + 1])*chivec.resid[i + 1]
            }
                ## Probability of an individual being undetected.
            log.p.unseen <- log(sum(((1 - pvec)*ptrvec.full + (1 - pvec)*chivec.resid*(1 - ptrvec.full))*pentvec))
            ## Likelihood contribution from the unseen individuals.
            nllike <- -lgamma(Ns + 1) + lgamma(Ns - nhist + 1) - (Ns - nhist) * log.p.unseen
            ## Constant so that we get the same log-likelihood as Robin's code.
            nllike <- nllike + lfactorial(nhist)
            ## Likelihood contributions from animals that were seen.
            for (f in 1:k){
                if (any(first.obs == f)){
                    ## (f for first-obs)
                    ## Extract the data records with first sighting on occasion f:
                    dat.f <- det.dat[first.obs==f, ,drop=F]
                    last.dat.f <- last.obs[first.obs==f]
                    n.f <- length(last.dat.f)
                    ## Calculate the psi parameters for this value of first-obs.
                    ## For a particular f, psif.vec is a vector of length f, where
                    ## psif.vec[i] = P(not seen on occasions i, ..., f-1 and survives to occasion f | alive at i)
                    psif.vec <- numeric(f)
                    psif.vec[f] <- 1
                    if(f > 1){
                        for( i in (f-1):1) {
                            psif.vec[i] = (1 - pvec[i]) * phivec[i] * psif.vec[i+1]
                        }
                    }
                    
                    ## All residents with first observation at time f have the same probability up to time f-1:
                    ## prob.to.f.m.1 = (pent1 * psif.1 + ... + pent.f * psif.f). Note that the mixture component
                    ## due to transients is zero, because transients don't survive.
                    prob.to.f.m.1 <- sum(pentvec[1:f] * psif.vec)
                    log.prob.to.f.m.1 <- log(prob.to.f.m.1)
                    
                    ## We need p(capt | transient) p(transient) + p(capt | resident) p(resident).
                    for (hst in 1:n.f){
                        last.hst <- last.dat.f[hst]
                        ## datvec.hst is the string of 0s and 1s for the single history hst, from the occasion f of its
                        ## first sighting to the occasion L-1 immediately previous to its last sighting.
                        datvec.hst <- dat.f[hst, f:last.hst]
                        p.hst <- pvec[f:last.hst]
                        #browser()
                        phi.hst <- phivec[f:(last.hst-1)]
                        ## For residents, joint probability of entry and obtaining all leading zeroes.

                        ## same as the prob.to.f.m.1
                        prob.resid.et.start <- pentvec[1:f]*(1 - ptrvec.full[1:f])*psif.vec
                        log.prob.resid.hst <- sum(datvec.hst*log(p.hst) +
                                                  (1 - datvec.hst)*log(1 - p.hst))
                        ## Survival between first and last detections if we have more than one.
                        if (sum(datvec.hst) > 1){
                            log.prob.resid.hst <- log.prob.resid.hst + sum(log(phi.hst))
                        }
                        prob.resid.hst <- exp(log.prob.resid.hst)
                        prob.resid.end <- chivec.resid[last.hst]
                        prob.resid.et <- prob.resid.et.start*prob.resid.hst*prob.resid.end
                        ## Capture histories for transients are impossible
                        ## unless they have a single detection, in which
                        ## case it must enter in the right session and be
                        ## detected.
                        prob.trans.et <- numeric(f)
                        if (sum(dat.f[hst, ]) == 1){
                            prob.trans.et[f] <- pentvec[f]*ptrvec.full[f]*pvec[f]
                        }
                        ## Putting it all together for this animal.
                        nllike <- nllike - log(sum(prob.resid.et + prob.trans.et))
                    }
                }
            }
 
            } else if (liktype == "cpp"){
                nllike <- popan_ll(k, phivec, pvec, ptrvec.full, pentvec, Ns, nhist, first.obs, last.obs, det.dat)
            } 
            nll <- nllike
            out <- get(out)
            out
        }
        ## Overall negative-log-likelihood is the sum across groups:
        nll <- sum(sapply(1:ngp, onegroup.func, out = out))
        if(printit){
            #for (i in 1:length(pars)){
            #    cat(names(pars)[i], ": ", round(pars[i], 3), ", ", sep = "")
            #}
            cat("NLL:", nll, "\n")
        }
        if (out == "nll"){
            out <- nll
        } else {
            out <- sapply(1:ngp, onegroup.func, out = out)
        }
        out
    }


    ## --------------------------------------------------------------------------------------------------
    ## Fit popan model and return:
    ## -----------------------------------------------------------------------------------------------------
    fit <- nlminb(
        start = startvec,
        objective = negloglik.func,
        lower = lowervec,
        upper = uppervec,
        control=list(iter.max=1e3, eval.max=1e3, rel.tol = 1e-15))
    pents <- negloglik.func(pars = fit$par, out = "pentvec")
    rownames(pents) <- NULL
    phis <- negloglik.func(pars = fit$par, out = "phivec")
    rhos <- negloglik.func(pars = fit$par, out = "rhovec")
    ps <- negloglik.func(pars = fit$par, out = "pvec")
    Ns <- negloglik.func(pars = fit$par, out = "Ns")
    ptrs <- negloglik.func(pars = fit$par, out = "ptrvec")
    ENs <- matrix(0, nrow = k, ncol = ngp)
    ENs[1, ] <- Ns*pents[1, ]
    ENs[2, ] <- (phis[1, ] + rhos[1, ])*(1 - ptrs[1, ])*ENs[1, ]
    for (i in 3:k){
        ENs[i, ] <- (phis[i - 1, ] + rhos[i - 1, ])*(phis[i - 2, ] + (1 - ptrs[i - 1, ])*rhos[i - 2, ])*ENs[i - 1,]/(phis[i - 2, ] + rhos[i - 2, ])
    }
    out <- list(fit = fit, Ns = Ns, phis = phis, rhos = rhos, ps = ps, ptrs = ptrs, pents = pents, ENs = ENs)
    class(out) <- "popan"
    out
}
