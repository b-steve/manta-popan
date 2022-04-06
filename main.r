## Function to fit a POPAN model.

## Arguments:
## captlist: A list with components for each group. Each component
##          contains capture histories.
## model.list: Model specifications for recruitment (b), survival
##             (phi), and detection probability (p). The word
##             "occasion" is reserved to indicate a separately
##             estimate for each occasion is required. For example,
##             model.list = list(b = ~ 1, phi = ~ x, p = ~ occasion)
##             means that recruitment is assumed to be constant over
##             time, survival depends on the covariate x, and
##             detection probability is separately estimated for each
##             occasion.
## group.pars: A named list, where names are parameter names (b, phi,
##             or p), and components are logical. If a component is
##             TRUE, then that parameter is shared between groups. For
##             example group.pars = list(b = FALSE, phi = FALSE, p =
##             TRUE) means that recruitment and survival are
##             separately estimated for each group, but detection
##             probabilities are assumed to be the same between
##             groups.
## group.effect: A named list, where names are parameter names (b,
##               phi, or p), and components are logical. If a
##               component is TRUE, then a group effect is fitted. If
##               FALSE, then parameters are separately estimated
##               between groups. For example, say model.list includes
##               the component p = ~ x, and group.pars includes
##               component p = FALSE. If group.effect includes
##               component p = TRUE, then the coefficient for
##               covariate x is used for all groups, but group is also
##               included as a categorical explanatory variable to
##               determine detection probabilities. If FALSE, then
##               the coefficent for x is separately estimated for each
##               group.
## random.start: Either TRUE or FALSE. If TRUE, then random start
##               values are used.
## startvec: A vector of specified start values (not recommended
##           because figuring out the model parameters is difficult).
## n.attempts: Number of attempts to fit the model using randomly
##             selected start values.
## n.cores: Number of cores for parallel processing when n.cores is
##          greater than 1.
fit.popan <- function(captlist, model.list = NULL, group.pars = NULL, group.effect = NULL, df = NULL, printit = FALSE, random.start = FALSE, startvec = NULL, n.attempts = 1, n.cores = 1){
    ## Saving function arguments.
    args <- list(captlist = captlist, model.list = model.list,
                 group.pars = group.pars, group.effect = group.effect,
                 df = df, printit = printit, random.start = random.start,
                 startvec = startvec, n.attempts = n.attempts, n.cores = 1)
    ## Detecting whether or not we are fitting a transience model.
    transience <- !is.null(model.list$ptr)
    if (n.attempts == 1){
        random.scale <- 1
        ## If captlist is a matrix, turn it into a list for consistency.
        if (!is.list(captlist)){
            captlist <- list(captlist)
        }
        ## Number of different groups.
        n.groups <- length(captlist)
        ## Number of occasions (needs to be the same across groups).
        k  <- ncol(captlist[[1]])
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
        ## Logical value for whether or not to put a group effect on b.
        b.effect <- group.effect[["b"]]
        if (is.null(b.effect)){
            b.effect <- FALSE
        }
        ## Error if we have both effect and not grouping.
        if (!b.group & b.effect){
            stop("You can't fit a group effect on parameter b if effects are not grouped.")
        }
        ## Put an additive group effect.
        if (group.effect[["b"]]){
            n.b.par <- n.b.par + 1
            b.par.names[[2]][1] <- "b1.2"
        }    
        ## Start values for b parameters.
        b.startvec <- numeric(sum(sapply(b.par.names, length)))
        names(b.startvec) <- c(b.par.names, recursive = TRUE)
        if (random.start){
            b.startvec[(substr(names(b.startvec), 1, 3) == "b1.")] <-
                rnorm(sum(substr(names(b.startvec), 1, 3) == "b1."), log(0.25), 0.2*random.scale)
            b.startvec[(substr(names(b.startvec), 1, 3) != "b1.")] <-
                rnorm(sum(substr(names(b.startvec), 1, 3) != "b1."), 0, 0.5*random.scale)
        }
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
        ## Logical value for whether or not to put a group effect on phi.
        phi.effect <- group.effect[["phi"]]
        if (is.null(phi.effect)){
            phi.effect <- FALSE
        }
        ## Error if we have both effect and not grouping.
        if (!phi.group & phi.effect){
            stop("You can't fit a group effect on parameter phi if effects are not grouped.")
        }
        ## Put an additive group effect.
        if (phi.effect){
            n.phi.par <- n.phi.par + 1
            phi.par.names[[2]][1] <- "phi1.2"
        }
        ## Start values for phi parameters.
        phi.startvec <- numeric(sum(sapply(phi.par.names, length)))
        names(phi.startvec) <- c(phi.par.names, recursive = TRUE)
        
        if (random.start){
            phi.startvec[(substr(names(phi.startvec), 1, 5) == "phi1.")] <-
                rnorm(sum(substr(names(phi.startvec), 1, 5) == "phi1."), qlogis(0.8), 0.5*random.scale)
            phi.startvec[(substr(names(phi.startvec), 1, 5) != "phi1.")] <-
                rnorm(sum(substr(names(phi.startvec), 1, 5) != "phi1."), 0, 0.25*random.scale)
        }
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
        ## Logical value for whether or not to put a group effect on p.
        p.effect <- group.effect[["p"]]
        if (is.null(p.effect)){
            p.effect <- FALSE
        }
        ## Error if we have both effect and not grouping.
        if (!p.group & p.effect){
            stop("You can't fit a group effect on parameter p if effects are not grouped.")
        }
        ## Put an additive group effect.
        if (group.effect[["p"]]){
            n.p.par <- n.p.par + 1
            p.par.names[[2]][1] <- "p1.2"
        }
        ## Start values for p parameters.
        p.startvec <- numeric(sum(sapply(p.par.names, length)))
        names(p.startvec) <- c(p.par.names, recursive = TRUE)
        p.startvec[(substr(names(p.startvec), 1, 3) == "p1.")] <- qlogis(0.1)
        if (random.start){
            p.startvec[(substr(names(p.startvec), 1, 3) == "p1.")] <-
                rnorm(sum(substr(names(p.startvec), 1, 3) == "p1."), qlogis(0.1), 0.5*random.scale)
            p.startvec[(substr(names(p.startvec), 1, 3) != "p1.")] <-
                rnorm(sum(substr(names(p.startvec), 1, 3) != "p1."), 0, 2)
        }








        
        ## Doing all the same stuff with ptr.
        ## Model for survival parameters.
        ptr.model <- model.list[["ptr"]]
        if (is.null(ptr.model)){
            ptr.model <- ~ occasion
        }
        ## Creating a function to model transience probability parameters. Note that it
        ## doesn't use the covariates for the final occasion.
        ptr.obj <- cov.func(ptr.model, df[-k, , drop = FALSE], plogis)
        ## The function.
        ptr.func <- ptr.obj[[1]]
        ## Number of parameters.
        n.ptr.par <- ptr.obj[[2]]
        ## Setting defaults so Rachel's function knows how many parameters there are.
        formals(ptr.func)$pars <- rep(0, n.ptr.par)
        ## Parameter names for the groups.
        ptr.par.names <- vector(mode = "list", length = n.groups)
        ## Logical value for whether or not to share ptr parameters across groups.
        ptr.group <- group.pars[["ptr"]]
        if (is.null(ptr.group)){
            ptr.group <- TRUE
        }
        for (i in 1:n.groups){
            ptr.par.names[[i]] <- paste0("ptr", 1:n.ptr.par, ".", ifelse(ptr.group, 1, i))
        }
        ## Logical value for whether or not to put a group effect on ptr.
        ptr.effect <- group.effect[["ptr"]]
        if (is.null(ptr.effect)){
            ptr.effect <- FALSE
        }
        ## Error if we have both effect and not grouping.
        if (!ptr.group & ptr.effect){
            stop("You can't fit a group effect on parameter ptr if effects are not grouped.")
        }
        ## Put an additive group effect.
        if (ptr.effect){
            n.ptr.par <- n.ptr.par + 1
            ptr.par.names[[2]][1] <- "ptr1.2"
        }
        ## Start values for ptr parameters.
        ptr.startvec <- numeric(sum(sapply(ptr.par.names, length)))
        names(ptr.startvec) <- c(ptr.par.names, recursive = TRUE)
        ptr.startvec[(substr(names(ptr.startvec), 1, 5) == "ptr1.")] <- qlogis(0.1)
        if (random.start){
            ptr.startvec[(substr(names(ptr.startvec), 1, 5) == "ptr1.")] <-
                rnorm(sum(substr(names(ptr.startvec), 1, 5) == "ptr1."), qlogis(0.1), 0.5*random.scale)
            ptr.startvec[(substr(names(ptr.startvec), 1, 5) != "ptr1.")] <-
                rnorm(sum(substr(names(ptr.startvec), 1, 5) != "ptr1."), 0, 2)
        }
        ## Start values for the N parameters.
        if (random.start){
            Ns.startvec <- c(Ns.1 = runif(1, nrow(captlist[[1]]), 3*nrow(captlist[[1]])),
                             Ns.2 = runif(1, nrow(captlist[[2]]), 3*nrow(captlist[[2]])))
        } else {
            Ns.startvec <- c(Ns.1 = 2*nrow(captlist[[1]]), Ns.2 = 2*nrow(captlist[[2]]))
        }
        ## Creating model object.
        model <- list()
        for (i in 1:n.groups){
            model[[i]] <- c(paste0("Ns.", i), b.par.names[[i]], phi.par.names[[i]], p.par.names[[i]], ptr.par.names[[i]][transience])
        }
        names(model) <- paste0("gp", 1:n.groups)
        ## Putting together the start values.
        if (is.null(startvec) | random.start){
            startvec <- c(Ns.startvec, b.startvec, phi.startvec, p.startvec, ptr.startvec[transience])
        }
        ## Fitting the model.
        if (transience){
            out <- popanGeneral.covs.fit.func.transience(captlist, k = k, birthfunc = b.func, phifunc = phi.func,
                                                         pfunc = p.func, ptrfunc = ptr.func,
                                                         model = model, startvec = startvec,
                                                         printit = printit)
        } else {
            out <- popanGeneral.covs.fit.func(captlist, k = k, birthfunc = b.func, phifunc = phi.func,
                                              pfunc = p.func, model = model, startvec = startvec,
                                              printit = printit)
        }
        out$args <- args
    } else {
        all.out <- vector(mode = "list", length = n.attempts)
        single.args <- args
        single.args$random.start <- TRUE
        single.args$n.attempts <- 1
        FUN <- function(i, args){
            ## First attempt is with default start values.
            if (i == 1) args$random.start <- FALSE
            try(do.call("fit.popan", args), silent = TRUE)
        }
        cluster <- makeCluster(n.cores)
        clusterEvalQ(cluster, {
            source("main.r")
        })
        all.out <- parLapplyLB(cluster, 1:(n.attempts + 1), FUN, args = single.args)
        stopCluster(cluster)
        lls <- sapply(all.out, function(x) -x$fit$objective)
        out <- all.out[[which(lls == max(lls))[1]]]
        out$args <- args
        out$all.lls <- lls
    }
    out
}

## Function to simulate new data based on estimates from a previously fitted model.

## Arguments:
## fit: A model object returned by fit.popan() from which to simulate new data.
sim.popan <- function(fit){
    Nsvec <- fit$Ns
    ngp <- length(Nsvec)
    phi.mat <- fit$phis
    rho.mat <- fit$rhos
    p.mat <- fit$ps
    k <- nrow(p.mat)
    phiList <- vector(mode = "list", length = ngp)
    rhoList <- vector(mode = "list", length = ngp)
    pList <- vector(mode = "list", length = ngp)
    for (i in 1:ngp){
        phiList[[i]] <- phi.mat[, i]
        rhoList[[i]] <- rho.mat[, i]
        pList[[i]] <- p.mat[, i]
    }
    out <- popanGeneral.covs.sim.func(k, Nsvec, phiList, rhoList, pList)
}

## Function to calculate the (Q)AIC for a POPAN model fit.

## Arguments:
## fit: A model object returned by fit.popan().
## chat: An estimate of c-hat for QAIC.
AIC.popan <- function(object, chat = 1, ...){
    ## Log-likleihood.
    ll <- -object$fit$objective
    ## Number of paramters.
    k <- length(object$fit$par)
    ## AIC.
    -2*ll/chat + 2*k
}

## Function to carry out a parametric bootstrap.

## Arguments:
## fit: A model object returned by fit.popan() to bootstrap.
## n.boots: The number of bootstrap iterations.
## n.cores: The number of cores for parallel computing. 
boot.popan <- function(fit, n.boots = 100, n.cores = 1){
    boots <- vector(mode = "list", length = n.boots)
    args.single <- fit$args
    args.full <- vector(mode = "list", length = n.boots)
    for (i in 1:n.boots){
        args.full[[i]] <- args.single
        args.full[[i]]$captlist <- sim.popan(fit)
    }
    fit$boots <- par.fit.popan(n.cores, arg.list = args.full)
    fit
}

## Function to summarise a bootstrap model fit, providing point
## estimates, standard errors, and 95% CIs for any function of model
## parameters.

## Arguments:
## boot.fit: A model object returned by boot.popan().
## par.fun: A function that computes a function of parameters to
##          summarise, with the primary input being the model object.
summary.popan <- function(boot.fit, par.fun = function(fit) fit$fit$par){
    ests <- par.fun(boot.fit)
    boots.par <- sapply(boot.fit$boots, par.fun)
    ses <- apply(boots.par, 1, sd)
    lower.ci <- apply(boots.par, 1, quantile, probs = 0.025)
    upper.ci <- apply(boots.par, 1, quantile, probs = 0.975)
    out <- matrix(0, nrow = length(ests), ncol = 4)
    out[, 1] <- ests
    out[, 2] <- ses
    out[, 3] <- lower.ci
    out[, 4] <- upper.ci
    colnames(out) <- c("Estimate", "Std Error", "Lower CI", "Upper CI")
    rownames(out) <- names(ests)
    out
}

## Function to carry out model averaging via a nonparametric
## bootstrap.

## Arguments:

## fits: A list, where each component is a candidate model object
##       returned by fit.popan().
## n.boots: The number of bootstrap iterations.
## chat: The c-hat value to use for QAIC. The default is 1, which
##       corresponds to AIC.
## n.cores: The number of cores for parallel computing.
## progress.bar: Logical. If TRUE, then a progress bar appears in the console.
boot.ma.popan <- function(fits, n.boots = 10, chat = 1, n.cores = 1, progress.bar = TRUE){
    n.fits <- length(fits)
    aics <- sapply(fits, AIC, chat = chat)
    ## Getting point estimates using weighted AIC model averaging from
    ## the original fits.
    min.aic <- min(aics)
    w <- exp((min.aic - aics)/2)
    w <- w/sum(w)
    phis.weighted <- 0*fits[[1]]$phis
    rhos.weighted <- 0*fits[[1]]$rhos
    ps.weighted <- 0*fits[[1]]$ps
    pents.weighted <- 0*fits[[1]]$pents
    ENs.weighted <- 0*fits[[1]]$ENs
    for (j in 1:n.fits){
        phis.weighted <- phis.weighted + w[j]*fits[[j]]$phis
        rhos.weighted <- rhos.weighted + w[j]*fits[[j]]$rhos
        ps.weighted <- ps.weighted + w[j]*fits[[j]]$ps
        pents.weighted <- pents.weighted + w[j]*fits[[j]]$pents
        ENs.weighted <- ENs.weighted + w[j]*fits[[j]]$ENs
    }
    out.weighted.ests <- list(phis = phis.weighted,
                              rhos = rhos.weighted,
                              ps = ps.weighted,
                              pents = pents.weighted,
                              ENs = ENs.weighted,
                              ENs.tot = apply(ENs.weighted, 1, sum))
    ## Getting original capture histories.
    captlist <- fits[[1]]$args$captlist
    ## Initialising args list.
    args.boot <- lapply(fits, function(x) x$args)
    ## Initialising output list.
    out.best <- vector(mode = "list", length = n.boots)
    out.weighted <- vector(mode = "list", length = n.boots)
    out.ests.boot <- vector(mode = "list", length = n.boots)
    ## Setting up progress bar.
    if (progress.bar){
        pb <- txtProgressBar(min = 0, max = n.boots, style = 3)
    }
    for (i in 1:n.boots){
        captlist.boot <- vector(mode = "list", length = length(captlist))
        names(captlist.boot) <- names(captlist)
        for (j in 1:length(captlist)){
            ## Getting bootstrap sample of each group by sampling with
            ## replacement.
            captlist.boot[[j]] <- captlist[[j]][sample(nrow(captlist[[j]]), replace = TRUE), ]
        }
        ## Putting the bootstrap data into the components of the args list.
        for (j in 1:n.fits){
            args.boot[[j]]$captlist <- captlist.boot
            args.boot[[j]]$n.attempts <- 1
            args.boot[[j]]$startvec <- fits[[j]]$fit$par
        }
        ## Fitting each of the models to the bootstrapped data set.
        fits.boot <- par.fit.popan(n.cores = n.cores, arg.list = args.boot)
        aics.boot <- sapply(fits.boot, AIC, chat = chat)
        ## Finding the best and saving results for the Buckland approach.
        which.model <- which(aics.boot == min(aics.boot))
        out.best[[i]] <- list(phis = fits.boot[[which.model]]$phis,
                              rhos = fits.boot[[which.model]]$rhos,
                              ps = fits.boot[[which.model]]$ps,
                              pents = fits.boot[[which.model]]$pents,
                              ENs = fits.boot[[which.model]]$ENs,
                              ENs.tot = apply(fits.boot[[which.model]]$ENs, 1, sum),
                              model.no = which.model)
        ## Alternatively, doing AIC weighting.
        min.aic.boot <- min(aics.boot)
        w.boot <- exp((min.aic.boot - aics.boot)/2)
        w.boot <- w.boot/sum(w.boot)
        phis.weighted.boot <- 0*fits.boot[[1]]$phis
        rhos.weighted.boot <- 0*fits.boot[[1]]$rhos
        ps.weighted.boot <- 0*fits.boot[[1]]$ps
        pents.weighted.boot <- 0*fits.boot[[1]]$pents
        ENs.weighted.boot <- 0*fits.boot[[1]]$ENs
        for (j in 1:n.fits){
            phis.weighted.boot <- phis.weighted.boot + w[j]*fits.boot[[j]]$phis
            rhos.weighted.boot <- rhos.weighted.boot + w[j]*fits.boot[[j]]$rhos
            ps.weighted.boot <- ps.weighted.boot + w[j]*fits.boot[[j]]$ps
            pents.weighted.boot <- pents.weighted.boot + w[j]*fits.boot[[j]]$pents
            ENs.weighted.boot <- ENs.weighted.boot + w[j]*fits.boot[[j]]$ENs
        }
        out.weighted[[i]] <- list(phis = phis.weighted.boot,
                                  rhos = rhos.weighted.boot,
                                  ps = ps.weighted.boot,
                                  pents = pents.weighted.boot,
                                  ENs = ENs.weighted.boot,
                                  ENs.tot = apply(ENs.weighted.boot, 1, sum))
        ## Saving all the fits for weighted bootstrap.
        out.ests.boot[[i]] <- lapply(fits.boot, function(x) x[c("phis", "rhos", "ps", "pents", "ENs")])
        ## Updating progress bar.
        if (progress.bar){
            setTxtProgressBar(pb, i)
        }
    }
    ## Closing progress bar.
    if (progress.bar){
        close(pb)
    }
    out <- list(best = out.best, weighted = out.weighted, weighted.ests = out.weighted.ests, ests.boot = out.ests.boot)
    class(out) <- "ma.popan"
    out
}

## Function to summarise model averaging results, providing point
## estimates, standard errors, and 95% CIs for any function of model
## parameters.

## Arguments:
## fit.ma: An object returned by boot.ma.popan().

## method: The method used to calculate summary output. If "weighted",
##         then on each bootstrap iteration the parameters retained
##         are those averaged across all candidate models. If "best",
##         then on each bootstrap iteration the parameters retained
##         are those from the model with the most AIC support.
## pars: The annual parameters to summarise. Options include "ENs",
##       "ps", "rhos", or "phis". Alternatively, the user may specify
##       their own function of mparameters to summarise using par.fun
## diffs: Logical. If TRUE, the parameters reported are all pairwise
##        differences between the requested parameters.
## par.fun: A function that computes a function of parameters to
##          summarise, with the primary input being the model object.
## par.fun.p: Logical. If TRUE, p-values are included to test the null
##            hypothesis that the true parameter value is equal to
##            zero. Defaults to FALSE because for many parameters it
##            is not meaningfull to test a null hypothesis of zero.
summary.ma.popan <- function(fit.ma, method = "best", pars = "ENs", groups = NULL, diffs = FALSE, par.fun = NULL, par.fun.p = FALSE){
    ## Bundling up all the arguments.
    args <- c(as.list(environment()))
    ## Total number of groups.
    n.groups <- ncol(fit.ma$best[[1]]$ENs)
    ## If group argument not provided, then do both groups.
    if (is.null(par.fun)){
        par.fun <- function(x, group) x[[pars]][, group]
        if (is.null(groups)){
            groups <- 1:n.groups
        }
    } else {
        groups <- NA
    }
    ## Recursive function because why not.
    if (length(groups) > 1){
        out <- vector(mode = "list", length = length(groups))
        names(out) <- paste0("group", groups)
        for (i in seq_along(groups)){
            args$groups <- groups[i]
            out[[i]] <- do.call("summary", args)
        }
    } else {
        boots.par <- sapply(fit.ma[[method]], par.fun, group = groups)
        if (!is.matrix(boots.par)){
            boots.par <- matrix(boots.par, nrow = 1)
        }
        if (diffs){
            boots.par <- apply(boots.par, 2, function(x) outer(x, x, `-`))
        }
        if (method == "weighted"){
            ests <- apply(boots.par, 1, mean)
        } else if (method == "best"){
            ests <- apply(boots.par, 1, mean)
        } else {
            stop("Argument 'method' must be \"best\" or \"weighted\"")
        }
        ses <- apply(boots.par, 1, sd)
        lower.ci <- apply(boots.par, 1, quantile, probs = 0.025)
        upper.ci <- apply(boots.par, 1, quantile, probs = 0.975)
        out <- matrix(0, nrow = length(ests), ncol = ifelse(diffs | par.fun.p, 5, 4))
        if (diffs | par.fun.p){
            ps <- apply(boots.par, 1, boot.p)
            out[, 5] <- ps
        }
        out[, 1] <- ests
        out[, 2] <- ses
        out[, 3] <- lower.ci
        out[, 4] <- upper.ci
        colnames(out) <- c("Estimate", "Std Error", "Lower CI", "Upper CI", "P-value"[diffs | par.fun.p])
    }
    out
}

## A wrapper function for the manta analyses, which carries out the
## following:

## - Fits a set of candidate models. Candidate specifications for each
##   demographic parameter (recruitment per capita, b; survival
##   probability, phi; and detection probability, p) are as follows:
##    - ~ 1
##      - constant for all years and across both sexes.
##    - ~ mei
##      - depends on annual MEI, but no differences between sexes.
##    - ~ occasion.
##      - separate estimates for each occasion, but no differences
##        between sexes.
##    - ~ sex
##      - constant for all years, but separately estimated for males
##        and females.
##    - ~ mei + sex
##      - depends on annual MEI and sex.
##    - ~ mei * sex
##      - depends on annual MEI, and relationship with MEI is
##        separately estimated for males and females.
##    - ~ occasion + sex
##      - separate estimates for each occasion, with an additive sex
##        effect so that the effect of sex is constant across all
##        years.
##    - ~ occasion * sex
##      - separate estimates for each occasion and sex.
##    - So in total there are 8 possible submodels for 3 different
##      demographic parameters, so the total number of candidate
##      models is 8^3 = 512
## - Carries out a model averaging procedure using a bootstrap.

## Arguments:
## captlist: A list with components for each group. Each component
##          contains capture histories.
## mei: The annual MEI values.
## n.boots: The number of bootstrap iterations.
## chat: The c-hat value to use for QAIC. The default is 1, which
##       corresponds to AIC.
## AIC.cutoff: A cutoff value to specify which models should be
##             considered when it comes to model averaging. Defaults
##             to 10, which means only models within 10 AIC units of
##             the 'best' model are considered. The more models that
##             are considered, the longer the function will take,
##             because all models under consideration are fitted for
##             each bootstrap iteration.
## random.start, n.atempts: Passed to fit.popan().
## n.cores: The number of cores for parallel computing.
## verbose: Logical. If TRUE, progress is printed to the console.

manta.ma.wrap <- function(captlist, mei, chat = 1, n.boots = 100, AIC.cutoff = 10, random.start = FALSE,
                          n.attempts = 1, n.cores = 1, verbose = TRUE){
    ## Turning mei into a data frame because that's what the fitting function needs.
    mei <- data.frame(mei)
    ## The object args is a list, where each component corresponds to a
    ## model to fit, and it contains the required arguments of fit.popan()
    ## for that model.
    args <- list()
    k <- 1
    for (model.b in 1:3){
        ## Model statement for recruitment.
        if (model.b == 1){
            group.b.set <- 2:3
            m.b <- ~ 1
        } else {
            group.b.set <- 1:3
            if (model.b == 2){
                m.b <- ~ mei
            } else if (model.b == 3){
                m.b  <- ~ occasion
            }
        }
        for (model.phi in 1:3){
            ## Model statement for survival.
            if (model.phi == 1){
                group.phi.set <- 2:3
                m.phi <- ~ 1
            } else {
                group.phi.set <- 1:3
                if (model.phi == 2){
                    m.phi <- ~ mei
                } else if (model.phi == 3){
                    m.phi  <- ~ occasion
                }
            }
            for (model.p in 1:3){
                ## Model statement for detection.
                if (model.p == 1){
                    group.p.set <- 2:3
                    m.p <- ~ 1
                } else {
                    group.p.set <- 1:3
                    if (model.p == 2){
                        m.p <- ~ mei
                    } else if (model.p == 3){
                        m.p  <- ~ occasion
                    }
                }
                for (group.b in group.b.set){
                    ## Parameter grouping for recruitment.
                    if (group.b == 1){
                        gp.b <- FALSE
                        ge.b <- FALSE
                    } else if (group.b == 2){
                        gp.b <- TRUE
                        ge.b <- TRUE
                    } else if (group.b == 3){
                        gp.b <- TRUE
                        ge.b <- FALSE
                    }
                    for (group.phi in group.phi.set){
                        ## Parameter grouping for survival.
                        if (group.phi == 1){
                            gp.phi <- FALSE
                            ge.phi <- FALSE
                        } else if (group.phi == 2){
                            gp.phi <- TRUE
                            ge.phi <- TRUE
                        } else if (group.phi == 3){
                            gp.phi <- TRUE
                            ge.phi <- FALSE
                        }
                        for (group.p in group.p.set){
                            ## Parameter grouping for detection.
                            if (group.p == 1){
                                gp.p <- FALSE
                                ge.p <- FALSE
                            } else if (group.p == 2){
                                gp.p <- TRUE
                                ge.p <- TRUE
                            } else if (group.p == 3){
                                gp.p <- TRUE
                                ge.p <- FALSE
                            }
                            args[[k]] <- list(captlist,
                                              model.list = list(b = m.b,
                                                                phi = m.phi,
                                                                p = m.p),
                                              group.pars = list(b = gp.b,
                                                                phi = gp.phi,
                                                                p = gp.p),
                                              group.effect = list(b = ge.b,
                                                                  phi = ge.phi,
                                                                  p = ge.p),
                                              df = mei, random.start = random.start,
                                              n.attempts = n.attempts, n.cores = 1)
                            k <- k + 1        
                        }
                    }
                }
            }
        }
    }
    ## Fitting all the models.
    if (verbose){
        message("Fitting models...")
    }
    fits <- par.fit.popan(n.cores, arg.list = args)
    ## Calculating AICs for all the models.
    fits.AICs <- sapply(fits, AIC, chat = chat)
    ## Keeping only those within the AIC cutoff.
    best.fits <- fits[fits.AICs - min(fits.AICs) <= AIC.cutoff]
    if (n.boots == 0){
        fit.ma <- NULL
    } else {
        ## Carrying out model averaging.
        if (verbose){
            message("Bootstrapping...")
        }
        fit.ma <- boot.ma.popan(best.fits, n.boots = n.boots, chat = chat,
                                n.cores = n.cores, progress.bar = verbose)
    }
    ## Returning object as output.
    list(best.fits = best.fits, fit.ma = fit.ma)
}

## Carrying out goodness-of-fit using R2ucare.
##
## Arguments:
## captlist: A list with components for each group. Each component
##          contains capture histories.
popan.gof <- function(captlist){
    n.groups <- length(captlist)
    ## Creating list to fill with test output.
    out <- vector(mode = "list", length = n.groups + 1)
    if (!is.null(names(captlist))){
        ## Using captlist group names if they are available.
        names(out) <- c(names(captlist), "combined")
    } else {
        ## Otherwise just call them "group1", "group2", etc.
        names(out) <- c(paste0("group", 1:n.groups), "combined")
    }
    for (i in 1:n.groups){
        ## Carry out the four tests for each group.
        out[[i]] <- vector(mode = "list", length = 6)
        names(out[[i]]) <- c("test2cl", "test2ct", "test3sr", "test3sm", "overall", "marray")
        captlist.freq <- rep(1, nrow(captlist[[i]]))
        out[[i]][[1]] <- test2cl(captlist[[i]], captlist.freq)
        out[[i]][[2]] <- test2ct(captlist[[i]], captlist.freq)
        out[[i]][[3]] <- test3sr(captlist[[i]], captlist.freq)
        out[[i]][[4]] <- test3sm(captlist[[i]], captlist.freq)
        out[[i]][[5]] <- overall_CJS(captlist[[i]], captlist.freq)
        out[[i]][[6]] <- marray(captlist[[i]], captlist.freq)
    }
    ## Combine the overall CJS test across all groups.
    overall.chi2 <- sum(sapply(out[1:n.groups], function(x) x[[5]]$chi2))
    overall.df <- sum(sapply(out[1:n.groups], function(x) x[[5]]$degree_of_freedom))
    overall.p <- 1 - pchisq(overall.chi2, overall.df)
    overall.chat <- overall.chi2/overall.df
    out[[n.groups + 1]] <- data.frame(chi2 = overall.chi2, degree_of_freedom = overall.df,
                                      p_value = overall.p, chat = overall.chat)
    out
}

popan.plot <- function(fit.ma, year.start = 2009, year.end = 2019) {

  ## Extracting parameter values from the model fit
  # ENs female + male
  EN.mf <- summary(fit.ma, par.fun = function(x, ...) x$ENs[, 1] + x$ENs[, 2])
  # ENs female
  EN.f <- summary(fit.ma, pars = "ENs", groups = 1)
  # ENs male
  EN.m <- summary(fit.ma, pars = "ENs", groups = 2)
  # ps female
  p.f <- summary(fit.ma, pars = "ps", groups = 1)
  # ps male
  p.m <- summary(fit.ma, pars = "ps", groups = 2)
  # phis female
  phi.f <- summary(fit.ma, pars = "phis", groups = 1)
  # phis male
  phi.m <- summary(fit.ma, pars = "phis", groups = 2)
  # rhos female
  rho.f <- summary(fit.ma, pars = "rhos", groups = 1)
  # rhos male
  rho.m <- summary(fit.ma, pars = "rhos", groups = 2)
  
  n.start <- (year.start-2009)+1
  n.end <- (year.end-2009)+1
  max.seq <- (n.end-n.start)+1
  
  ## PLOTTING ####
  par(mfrow=c(3,2)) # dividing plots into 6 panels (row, col)
  par(mar = c(3,3,2,1)) # adjust graph margin: below, left, top, right
  par(mgp = c(1.8,0.4,0))
  
  ## PANEL A = ENs for total population (female+male) and for females only ####
  # total population (female+male)
  plot(EN.mf[n.start:n.end, 1],
       col = "orangered2",
       ylim = c(0, max(EN.mf[, 4])),
       type = "b", 
       pch = 15,
       xlab = "Year",
       ylab = "EN",
       main = "Expected population size",
       cex.main = 1.7, # size of axis title
       xaxt="n",
       cex.axis = 1.2, # sze of axis label
       cex = 1.6, # size of symbols
       tck=-.015, # size of ticks in axis
       cex.lab = 1.6) # size of axis title
  axis(1, at = seq(1,max.seq, by = 1), cex.axis=1.2,
       labels = c(year.start:year.end), tck=-.015)
  ## Adding dotted line for lower and upper CI limits
  lines(EN.mf[n.start:n.end, 3], lty = "longdash", col = "orangered2")
  lines(EN.mf[n.start:n.end, 4], lty = "longdash", col = "orangered2")
  
  par(new = T)
  
  # plotting female
  plot(EN.f[n.start:n.end, 1], 
       col = "goldenrod1",
       ylim = c(0, max(EN.mf[, 4])),
       type = "b", 
       pch = 1,
       xlab = "",
       ylab = "",
       main = "",
       cex.main = 1.7,
       xaxt="n",
       yaxt="n",
       cex.axis = 1.2,
       cex = 1.6,
       tck=-.015,
       cex.lab = 1.6)
  ## Adding dotted line for lower and upper CI limits
  lines(EN.f[n.start:n.end, 3], lty = "dashed", col = "goldenrod1")
  lines(EN.f[n.start:n.end, 4], lty = "dashed", col = "goldenrod1")
  
  # create legend for Panel A
  legend("topleft", legend= c("Females + Males", "Females"), col = c('orangered2', 'goldenrod1'),
         cex=1.6, bg="transparent", pch = c(15,1),
         box.lty=0)
  legend("bottomright", text.font=2, legend= "A",
         cex=2.5, bg="transparent",
         box.lty=0)
  box()
  
  ## PANEL B = ENs for total population (female+male) and for females only ####
  # total population (female+male)
  plot(EN.mf[n.start:n.end, 1],
       col = "orangered2",
       ylim = c(0, max(EN.mf[, 4])),
       type = "b", 
       pch = 15,
       xlab = "Year",
       ylab = "EN",
       main = "Expected population size",
       cex.main = 1.7, # size of axis title
       xaxt="n",
       cex.axis = 1.2, # sze of axis label
       cex = 1.6,
       tck=-.015,
       cex.lab = 1.6)
  axis(1, at = seq(1,max.seq, by = 1), cex.axis=1.2,
       labels = c(year.start:year.end), tck=-.015)
  ## Adding dotted line for lower and upper CI limits
  lines(EN.mf[n.start:n.end, 3], lty = "longdash", col = "orangered2")
  lines(EN.mf[n.start:n.end, 4], lty = "longdash", col = "orangered2")
  
  par(new = T)
  
  # plotting male
  plot(EN.m[n.start:n.end, 1], 
       col = "turquoise3",
       ylim = c(0, max(EN.mf[, 4])),
       type = "b", 
       pch = 2,
       xlab = "",
       ylab = "",
       main = "",
       cex.main = 1.7,
       xaxt="n",
       yaxt="n",
       cex.axis = 1.2,
       cex = 1.6,
       tck=-.015,
       cex.lab = 1.6)
  ## Adding dotted line for lower and upper CI limits
  lines(EN.m[n.start:n.end, 3], lty = "dashed", col = "turquoise3")
  lines(EN.m[n.start:n.end, 4], lty = "dashed", col = "turquoise3")
  
  # create legend for Panel B
  legend("topleft", legend= c("Females + Males", "Males"), col = c('orangered2', 'turquoise3'),
         cex=1.6, bg="transparent", pch = c(15,2),
         box.lty=0)
  legend("bottomright", text.font=2, legend= "B",
         cex=2.5, bg="transparent",
         box.lty=0)
  box()
  
  ## PANEL C = sighting probability (p) for female ####
  plot(p.f[n.start:n.end, 1],
       col = "goldenrod1",
       ylim = c(0, max(p.f[, 4])),
       type = "b", 
       pch = 1,
       xlab = "Year",
       ylab = "p",
       main = "Sighting probabilities",
       cex.main = 1.7, # size of axis title
       xaxt="n",
       cex.axis = 1.2, # sze of axis label
       cex = 1.6,
       tck=-.015,
       cex.lab = 1.6)
  axis(1, at = seq(1,max.seq, by = 1), cex.axis=1.2,
       labels = c(year.start:year.end), tck=-.015)
  ## Adding dotted line for lower and upper CI limits
  lines(p.f[n.start:n.end, 3], lty = "longdash", col = "goldenrod1")
  lines(p.f[n.start:n.end, 4], lty = "longdash", col = "goldenrod1")
  
  # create legend for Panel C
  legend("topleft", legend= "Females", col = "goldenrod1",
         pch = 1, cex=1.6, bg="transparent",
         box.lty=0)
  legend("bottomright", text.font=2, legend= "C",
         cex=2.5, bg="transparent",
         box.lty=0)
  box()
  
  ## PANEL D = sighting probability (p) for male ####
  plot(p.m[n.start:n.end, 1],
       col = "turquoise3",
       ylim = c(0, max(p.f[, 4])),
       type = "b", 
       pch = 1,
       xlab = "Year",
       ylab = "p",
       main = "Sighting probabilities",
       cex.main = 1.7, # size of axis title
       xaxt="n",
       cex.axis = 1.2, # sze of axis label
       cex = 1.6,
       tck=-.015,
       cex.lab = 1.6)
  axis(1, at = seq(1,max.seq, by = 1), cex.axis=1.2,
       labels = c(year.start:year.end), tck=-.015)
  ## Adding dotted line for lower and upper CI limits
  lines(p.m[n.start:n.end, 3], lty = "longdash", col = "turquoise3")
  lines(p.m[n.start:n.end, 4], lty = "longdash", col = "turquoise3")
  
  # create legend for Panel C
  legend("topleft", legend= "Males", col = "turquoise3",
         pch = 1, cex=1.6, bg="transparent",
         box.lty=0)
  legend("bottomright", text.font=2, legend= "D",
         cex=2.5, bg="transparent",
         box.lty=0)
  box()
  
  ## PANEL E: survival probabilities (phis)
  # phis for females
  phi.f
  last.row <- rep(NA, ncol(phi.f))
  phi.f.full <- rbind(phi.f, last.row)
  ## Creating a line for point estimates. The y-axis goes from 0 to the
  ## highest upper CI limit.
  plot(phi.f.full[n.start:n.end, 1], 
       col = "goldenrod1",
       ylim = c(0, 1),
       type = "b", 
       pch = 1,
       xlab = "Year",
       ylab = "phi",
       main = "Survival probabilities",
       cex.main = 1.7,
       xaxt="n",
       cex.axis = 1.2,
       cex = 1.6,
       cex.lab = 1.6)
  axis(1, at = seq(1,max.seq, by = 1), cex.axis=1.2,
       labels = c(year.start:year.end), tck=-.015)
  
  ## Adding dotted line for lower and upper CI limits
  lines(phi.f.full[n.start:n.end, 3], lty = "dashed", col = "goldenrod1")
  lines(phi.f.full[n.start:n.end, 4], lty = "dashed", col = "goldenrod1")
  
  par(new = T)
  
  # phis for males
  phi.m.full <- rbind(phi.m, last.row)
  plot(phi.m.full[n.start:n.end, 1],
       col = "turquoise3",
       ylim = c(0, 1),
       type = "b", 
       pch = 2,
       xlab = "",
       ylab = "",
       xaxt="n",
       cex = 1.6,
       yaxt="n")
  lines(phi.m.full[n.start:n.end, 3], lty = "dotted", col = "turquoise3")
  lines(phi.m.full[n.start:n.end, 4], lty = "dotted", col = "turquoise3")
  
  # create legend for Panel E
  legend("bottomleft", legend= c("Females", "Males"), col = c("goldenrod1", "turquoise3"),
         pch = c(1,2), cex=1.6, bg="transparent",
         box.lty=0)
  legend("bottomright", text.font=2, legend= "E",
         cex=2.5, bg="transparent",
         box.lty=0)
  box()
  
  
  ## PANEL F: per capita recruitment rate (rho)
  # phis for females
  rho.f
  last.row <- rep(NA, ncol(rho.f))
  rho.f.full <- rbind(rho.f, last.row)
  ## Creating a line for point estimates. The y-axis goes from 0 to the
  ## highest upper CI limit.
  plot(rho.f.full[n.start:n.end, 1], 
       col = "goldenrod1",
       ylim = c(0, 1),
       type = "b", 
       pch = 1,
       xlab = "Year",
       ylab = "rho",
       main = "Per capita recruitment rates",
       cex.main = 1.7,
       xaxt="n",
       cex.axis = 1.2,
       cex = 1.6,
       cex.lab = 1.6)
  axis(1, at = seq(1,max.seq, by = 1), cex.axis=1.2,
       labels = c(year.start:year.end), tck=-.015)
  ## Adding dotted line for lower and upper CI limits
  lines(rho.f.full[n.start:n.end, 3], lty = "dashed", col = "goldenrod1")
  lines(rho.f.full[n.start:n.end, 4], lty = "dashed", col = "goldenrod1")
  
  par(new = T)
  
  # rhos for males
  rho.m.full <- rbind(rho.m, last.row)
  plot(rho.m.full[n.start:n.end, 1],
       col = "turquoise3",
       ylim = c(0, 1),
       type = "b", 
       pch = 2,
       xlab = "",
       ylab = "",
       xaxt="n",
       cex = 1.6,
       yaxt="n")
  lines(rho.m.full[n.start:n.end, 3], lty = "dotted", col = "turquoise3")
  lines(rho.m.full[n.start:n.end, 4], lty = "dotted", col = "turquoise3")
  
  # create legend for Panel F
  legend("topleft", legend= c("Females", "Males"), col = c("goldenrod1", "turquoise3"),
         pch = c(1,2), cex=1.6, bg="transparent",
         box.lty=0)
  legend("bottomright", text.font=2, legend= "F",
         cex=2.5, bg="transparent",
         box.lty=0)
  box()
  
}


## Some packages.
library(mgcv)
library(parallel)
library(R2ucare)
source("internals.r")
