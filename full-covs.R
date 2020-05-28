popanGeneral.covs.fit.func <- function(dat, k=ncol(dat[[1]]), birthfunc = immigrationElNino.func, phifunc, pfunc,
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
    nphipar <- length(eval(formals(phifunc)[[1]]))
    nppar <- length(eval(formals(pfunc)[[1]]))
    parnamesTemplate <- c("Ns", paste0("b", 1:nbpar), paste0("phi", 1:nphipar), paste0("p", 1:nppar))
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
    lowervec[grep("phi", parsEstvec)] <- -Inf
    lowervec[grep("p", parsEstvec)] <- -Inf
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
    ## uppervec is 1 for phi and p parameters, Inf for Ns, bpar, and phipar parameters:
    uppervec <- rep(1-1e-6, npar)
    uppervec[grep("Ns", parsEstvec)] <- Inf
    uppervec[grep("b", parsEstvec)] <- Inf
    uppervec[grep("phi", parsEstvec)] <- Inf
    uppervec[grep("p", parsEstvec)] <- Inf
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

            ## Use the birthfunc supplied to convert bpars into rhovec:
            rhovec <- birthfunc(bpars)
            ## Use phifunc supplied to convert phipars into phivec.
            phivec <- phifunc(phipars)
            ## Use phifunc supplied to convert phipars into phivec.
            pvec <- pfunc(ppars)
            
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
            nll <- -sum(
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
            out <- get(out)
            out
        }
        ## Overall negative-log-likelihood is the sum across groups:
        nll <- sum(sapply(1:ngp, onegroup.func, out = out))
        if(printit){
            for (i in 1:length(pars)){
                cat(names(pars)[i], ": ", round(pars[i], 3), ", ", sep = "")
            }
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
        control=list(iter.max=1000, eval.max=5000))
    pents <- negloglik.func(pars = fit$par, out = "pentvec")
    rownames(pents) <- NULL
    phis <- negloglik.func(pars = fit$par, out = "phivec")
    rhos <- negloglik.func(pars = fit$par, out = "rhovec")
    ps <- negloglik.func(pars = fit$par, out = "pvec")
    Ns <- negloglik.func(pars = fit$par, out = "Ns")
    ENs <- matrix(0, nrow = k, ncol = ngp)
    ENs[1, ] <- Ns*pents[1, ]
    for (i in 2:k){
        ENs[i, ] <- phis[i - 1, ]*ENs[i - 1] + Ns*pents[i, ]
    }
    out <- list(fit = fit, Ns = Ns, phis = phis, rhos = rhos, ps = ps, pents = pents, ENs = ENs)
    class(out) <- "popan"
    out
}

## Function that generates a parameter modelling function (e.g., birthfunc above).
## - formula is a standard model formula.
## - df is a data frame of covariates.
## - invlink is the inverse-link function to apply.
cov.func <- function(formula, df, invlink = identity){
    ## Number of observations.
    n.obs <- nrow(df)
    ## Creating a fake response to trick gam() into creating our model matrix.
    resp <- rep(0, n.obs)
    ## Sticking the response on the LHS of the formula.
    full.formula <- as.formula(paste("resp", paste(as.character(formula), collapse = "")))
    ## Creating an object using gam() without fitting a model.
    fgam <- gam(full.formula, data = df, fit = FALSE)
    mm <- fgam$X
    colnames(mm) <- fgam$term.names
    n.par <- ncol(mm)
    def <- rep(0, n.par)
    ## Returning a function that creates derived parameters from coefficients.
    list(fun = 
             function(pars){
                 c(invlink(mm %*% pars))
             },
         n.par = n.par)
}

fit.popan <- function(caplist, model.list = NULL, group.pars = NULL, df = NULL, printit = FALSE){
    ## If caplist is a matrix, turn it into a list for consistency.
    if (!is.list(caplist)){
        caplist <- list(caplist)
    }
    ## Number of different groups.
    n.groups <- length(caplist)
    ## Number of occasions (needs to be the same across groups).
    k  <- ncol(caplist[[1]])
    # Setting up data frame.
    if (is.null(df)){
        df <- data.frame(occasion = as.factor(1:k))
    } else {
        df$occasion <- as.factor(1:k)
    }
    ## Model for birth parameters.
    b.model <- model.list[["b"]]
    if (is.null(b.model)){
        b.model <- ~ occasion
    }
    ## Creating a function to model birth parameters. Note that it
    ## doesn't use the covariates for the final occasion.
    b.obj <- cov.func(b.model, df[-k, , drop = FALSE], exp)
    ## The function.
    b.func <- b.obj[[1]]
    ## Number of parameters.
    n.b.par <- b.obj[[2]]
    ## Setting defaults so Rachel's function knows how many parameters there are.
    formals(b.func)$pars <- rep(0, n.b.par)
    ## Parameter names for the groups.
    b.par.names <- vector(mode = "list", length = n.groups)
    ## Logical value for whether or not to share b parameters across groups.
    b.group <- group.pars[["b"]]
    if (is.null(b.group)){
        b.group <- TRUE
    }
    for (i in 1:n.groups){
        b.par.names[[i]] <- paste0("b", 1:n.b.par, ".", ifelse(b.group, 1, i))
    }
    ## Start values for b parameters.
    b.startvec <- numeric(sum(sapply(b.par.names, length)))
    names(b.startvec) <- c(b.par.names, recursive = TRUE)
    ## Doing all the same stuff with phi.
    ## Model for survival parameters.
    phi.model <- model.list[["phi"]]
    if (is.null(phi.model)){
        phi.model <- ~ occasion
    }
    ## Creating a function to model survival parameters. Note that it
    ## doesn't use the covariates for the final occasion.
    phi.obj <- cov.func(phi.model, df[-k, , drop = FALSE], plogis)
    ## The function.
    phi.func <- phi.obj[[1]]
    ## Number of parameters.
    n.phi.par <- phi.obj[[2]]
    ## Setting defaults so Rachel's function knows how many parameters there are.
    formals(phi.func)$pars <- rep(0, n.phi.par)
    ## Parameter names for the groups.
    phi.par.names <- vector(mode = "list", length = n.groups)
    ## Logical value for whether or not to share phi parameters across groups.
    phi.group <- group.pars[["phi"]]
    if (is.null(phi.group)){
        phi.group <- TRUE
    }
    for (i in 1:n.groups){
        phi.par.names[[i]] <- paste0("phi", 1:n.phi.par, ".", ifelse(phi.group, 1, i))
    }
    ## Start values for phi parameters.
    phi.startvec <- numeric(sum(sapply(phi.par.names, length)))
    names(phi.startvec) <- c(phi.par.names, recursive = TRUE)

    ## Doing all the same stuff with p.
    ## Model for survival parameters.
    p.model <- model.list[["p"]]
    if (is.null(p.model)){
        p.model <- ~ occasion
    }
    ## Creating a function to model survival parameters.
    p.obj <- cov.func(p.model, df, plogis)
    ## The function.
    p.func <- p.obj[[1]]
    ## Number of parameters.
    n.p.par <- p.obj[[2]]
    ## Setting defaults so Rachel's function knows how many parameters there are.
    formals(p.func)$pars <- rep(0, n.p.par)
    ## Parameter names for the groups.
    p.par.names <- vector(mode = "list", length = n.groups)
    ## Logical value for whether or not to share p parameters across groups.
    p.group <- group.pars[["p"]]
    if (is.null(p.group)){
        p.group <- TRUE
    }
    for (i in 1:n.groups){
        p.par.names[[i]] <- paste0("p", 1:n.p.par, ".", ifelse(p.group, 1, i))
    }
    ## Start values for p parameters.
    p.startvec <- numeric(sum(sapply(p.par.names, length)))
    names(p.startvec) <- c(p.par.names, recursive = TRUE)
    p.startvec[(substr(names(p.startvec), 1, 2) == "p1")] <- qlogis(0.1)
    ## Start values for the N parameters.
    Ns.startvec <- c(Ns.1 = 1000, Ns.2 = 1200)
    ## Creating model object.
    model <- list(gp1 = c("Ns.1", b.par.names[[1]], phi.par.names[[1]], p.par.names[[1]]),
                  gp2 = c("Ns.2", b.par.names[[2]], phi.par.names[[2]], p.par.names[[2]]))
    ## Putting together the start values.
    startvec <- c(Ns.startvec, b.startvec, phi.startvec, p.startvec)
    ## Fitting the model.
    popanGeneral.covs.fit.func(caplist, k = k, birthfunc = b.func, phifunc = phi.func, pfunc = p.func,
                               model = model, startvec = startvec, printit = printit)
}

par.fit.popan <- function(n.cores, ..., arg.list = NULL){
    if (is.null(arg.list)){
        arg.list <- list(...)
    }
    n.fits <- length(arg.list)
    FUN <- function(i, arg.list){
        out <- try(do.call(fit.popan, arg.list[[i]]), silent = TRUE)
    }
    cluster <- makeCluster(n.cores)
    clusterEvalQ(cluster, {
        library(mgcv)
        source("Rfunc.R")
        source("full-covs.R") 
    })
    out <- parLapplyLB(cluster, 1:n.fits, FUN, arg.list = arg.list)
    stopCluster(cluster)
    out
}
