## Function to fit a POPAN model.

## Arguments:
## caplist: A list with components for each group. Each component
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
fit.popan <- function(caplist, model.list = NULL, group.pars = NULL, group.effect = NULL, df = NULL, printit = FALSE){
    ## Saving function arguments.
    args <- list(caplist = caplist, model.list = model.list,
                 group.pars = group.pars, group.effect = group.effect,
                 df = df, printit = printit)
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

    ###
    #n.b.par <- 2
    
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
    p.startvec[(substr(names(p.startvec), 1, 2) == "p1")] <- qlogis(0.1)
    ## Start values for the N parameters.
    Ns.startvec <- c(Ns.1 = 1000, Ns.2 = 1200)
    ## Creating model object.
    model <- list(gp1 = c("Ns.1", b.par.names[[1]], phi.par.names[[1]], p.par.names[[1]]),
                  gp2 = c("Ns.2", b.par.names[[2]], phi.par.names[[2]], p.par.names[[2]]))
    ## Putting together the start values.
    startvec <- c(Ns.startvec, b.startvec, phi.startvec, p.startvec)
    ## Fitting the model.
    out <- popanGeneral.covs.fit.func(caplist, k = k, birthfunc = b.func, phifunc = phi.func,
                                      pfunc = p.func, model = model, startvec = startvec,
                                      printit = printit)
    out$args <- args
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
        args.full[[i]]$caplist <- sim.popan(fit)
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
                              petns = pents.weighted,
                              ENs = ENs.weighted,
                              ENs.tot = apply(ENs.weighted, 1, sum))
    ## Getting original capture histories.
    caplist <- fits[[1]]$args$caplist
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
        caplist.boot <- vector(mode = "list", length = length(caplist))
        names(caplist.boot) <- names(caplist)
        for (j in 1:length(caplist)){
            ## Getting bootstrap sample of each group by sampling with
            ## replacement.
            caplist.boot[[j]] <- caplist[[j]][sample(nrow(caplist[[j]]), replace = TRUE), ]
        }
        ## Putting the bootstrap data into the components of the args list.
        for (j in 1:n.fits){
            args.boot[[j]]$caplist <- caplist.boot
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
                                  petns = pents.weighted.boot,
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
## diffs: Logical. If TRUE, the parameters reported are differences.
## par.fun: A function that computes a function of parameters to
##          summarise, with the primary input being the model object.
## par.fun.p: Logical. If TRUE, p-values are included to test the null
##            hypothesis that the true parameter value is equal to
##            zero. Defaults to FALSE because for many parameters it
##            is not meaningfull to test a null hypothesis of zero.
summary.ma.popan <- function(fit.ma, method = "weighted", pars = "ENs", diffs = FALSE, par.fun = NULL, par.fun.p = FALSE){
    if (is.null(par.fun)){
        par.fun <- function(x) x[[pars]]
    }
    boots.par <- sapply(fit.ma[[method]], par.fun)
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
## caplist: A list with components for each group. Each component
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
## n.cores: The number of cores for parallel computing.
## verbose: Logical. If TRUE, progress is printed to the console.

manta.ma.wrap <- function(captlist, mei, chat = 1, n.boots = 100, AIC.cutoff = 10, n.cores = 1, verbose = TRUE){
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
                                              df = covs)
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
    best.fits <- fits[fits.AICs - min(fits.AICs) <= 10]
    ## Carrying out model averaging.
    if (verbose){
        message("Bootstrapping...")
    }
    fit.ma <- boot.ma.popan(best.fits, n.boots = n.boots, chat = chat,
                            n.cores = n.cores, progress.bar = verbose)
    ## Returning object as output.
    fit.ma
}

## Some packages.
library(mgcv)
library(parallel)
library(R2ucare)
source("internals.r")
