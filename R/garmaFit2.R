garmaFit2 <-
function (formula = formula(data), order = c(0, 0), weights = NULL, 
    data = sys.parent(), family = NO(), alpha = 0.1, phi.start = NULL, 
    theta.start = NULL, tail = max(order), control = list()) 
{
    rqres <- function(pfun = "pNO", type = c("Continuous", "Discrete", 
        "Mixed"), censored = NULL, ymin = NULL, mass.p = NULL, 
        prob.mp = NULL, y = y, ...) {
    }
    body(rqres) <- eval(quote(body(rqres)), envir = getNamespace("gamlss"))
    hessian <- function(par) {
        npar <- length(par)
        epsilon <- 1e-04 * par
        Hessian = matrix(0, ncol = npar, nrow = npar)
        for (i in 1:npar) {
            for (j in 1:npar) {
                x1. <- x2. <- x3. <- x4. <- par
                x1.[i] <- x1.[i] + epsilon[i]
                x1.[j] <- x1.[j] + epsilon[j]
                x2.[i] <- x2.[i] + epsilon[i]
                x2.[j] <- x2.[j] - epsilon[j]
                x3.[i] <- x3.[i] - epsilon[i]
                x3.[j] <- x3.[j] + epsilon[j]
                x4.[i] <- x4.[i] - epsilon[i]
                x4.[j] <- x4.[j] - epsilon[j]
                Hessian[i, j] <- (garmaLL(x1.) - garmaLL(x2.) - 
                  garmaLL(x3.) + garmaLL(x4.))/(4 * epsilon[i] * 
                  epsilon[j])
            }
        }
        Hessian
    }
    createLags <- function(y, lag = 1, omit.na = FALSE) {
        lag1 <- function(y) {
            l <- c(NA, y[1:(length(y) - 1)])
            l
        }
        d <- matrix(0, nrow = length(y), ncol = lag)
        d[, 1] <- lag1(y)
        yname <- deparse(substitute(y))
        cname <- paste("lag1", yname, sep = "")
        if (lag == 1) {
            names(d[, 1]) <- c("lag1")
            if (omit.na) 
                d <- na.omit(d)
            return(d)
        }
        for (i in 2:lag) {
            d[, i] <- lag1(d[, i - 1])
            cname <- c(cname, paste("lag", paste(yname, i, sep = ""), 
                sep = ""))
        }
        colnames(d) <- cname
        if (omit.na) 
            d <- na.omit(d)
        d
    }
    garmaLL <- function(parm, save = FALSE) {
        beta <- parm[1:l.beta]
        Xb <- X %*% beta
        switch(case, {
            start <- l.beta + 1
            finish <- start + (order[1] - 1)
            phi <- parm[start:finish]
            finish <- finish + 1
            g2ylag <- createLags(g2y, order[1])
            Xblag <- createLags(Xb, order[1])
            g2y_Xb <- g2ylag - Xblag
            eta <- Xb + g2y_Xb %*% phi
            mu <- ifelse(is.na(eta), alpha, fam$mu.linkinv(eta))
        }, {
            start <- l.beta + 1
            finish <- start + (order[2] - 1)
            theta <- parm[start:finish]
            finish <- finish + 1
            g2ylag <- createLags(g2y, order[2])
            some <- Xb + g2ylag %*% theta
            some <- ifelse(is.na(some), 0, some)
            eta <- filter(some, -theta, "r", init = rep(0, order[2]))
            mu <- ifelse(is.na(eta), alpha, fam$mu.linkinv(eta))
        }, {
            start <- l.beta + 1
            finish1 <- start + (order[1] - 1)
            finish2 <- finish1 + 1 + (order[2] - 1)
            phi <- parm[start:finish1]
            theta <- parm[(finish1 + 1):finish2]
            finish <- finish2 + 1
            g2ylagar <- createLags(g2y, order[1])
            g2ylagma <- createLags(g2y, order[2])
            Xblagar <- createLags(Xb, order[1])
            Xblagma <- createLags(Xb, order[2])
            g2y_Xb <- g2ylagar - Xblagar
            some <- Xb + g2y_Xb %*% phi + g2ylagma %*% theta
            some <- ifelse(is.na(some), 0, some)
            eta <- filter(some, -theta, "r", init = rep(0, order[2]))
            mu <- ifelse(is.na(eta), alpha, fam$mu.linkinv(eta))
        })
        switch(nopar, {
            llh <- if (BItrue) -sum(w * PDF(y, mu = mu, bd = bd, 
                log = TRUE), na.rm = TRUE) else -sum(w * PDF(y, 
                mu = mu, log = TRUE), na.rm = TRUE)
        }, {
            sigma <- parm[finish]
            llh <- -sum(w * PDF(y, mu = mu, sigma = sigma, log = TRUE), 
                na.rm = TRUE)
        }, {
            sigma <- parm[finish]
            nu <- parm[finish + 1]
            llh <- -sum(w * PDF(y, mu = mu, sigma = sigma, nu = nu, 
                log = TRUE), na.rm = TRUE)
        }, {
            sigma <- parm[finish]
            nu <- parm[finish + 1]
            tau <- parm[finish + 2]
            llh <- -sum(w * PDF(y, mu = mu, sigma = sigma, nu = nu, 
                tau = tau, log = TRUE), na.rm = TRUE)
        })
        if (save) 
            return(list(lik = llh, mu = mu))
        llh
    }
    garmacall <- match.call()
    if (order[1] <= 0 && order[2] <= 0) 
        case <- 0
    if (order[1] > 0 && order[2] <= 0) 
        case <- 1
    if (order[1] <= 0 && order[2] > 0) 
        case <- 2
    if (order[1] > 0 && order[2] > 0) 
        case <- 3
    m0 <- gamlss(formula, family = family, data = data, trace = FALSE)
    cat("deviance of linear model= ", deviance(m0), "\n")
    if (case == 0) 
        return(m0)
    beta <- coef(m0)
    l.beta <- length(beta)
    y <- m0$y
    X <- m0$mu.x
    fam <- as.gamlss.family(family)
    nopar <- fam$nopar
    N <- length(y)
    if (any(fam$family %in% gamlss::.gamlss.bi.list)) {
        BItrue <- TRUE
        bd <- m0$bd
    } else (BItrue <- FALSE)
    fname <- fam$family[[1]]
    dfun <- paste("d", fname, sep = "")
    pfun <- paste("p", fname, sep = "")
    PDF <- eval(parse(text = dfun))
    CDF <- eval(parse(text = pfun))
    if (fam$mu.link == "identity") 
        ystar <- y
    if (fam$mu.link == "log") 
        ystar <- pmax(y, alpha)
    if (fam$mu.link == "logit") 
        ystar <- pmin(pmax(y, alpha), bd - alpha)/bd
    if (fam$mu.link == "inverse") 
        ystar <- y
    g2y <- fam$mu.linkfun(ystar)
    tailoff <- max(order)
    mtail <- max(tailoff, tail)
    w <- if (is.null(weights)) 
        rep(1, N)
    else weights
    if (any(w < 0)) 
        stop("negative weights not allowed")
    w[1:mtail] <- 0
    if ("mu" %in% names(fam$parameters)) {
        params <- c(beta = beta)
        lowerBounds <- c(beta = rep(-Inf, l.beta))
        upperBounds <- c(beta = rep(Inf, l.beta))
        switch(case, {
            if (length(phi.start) != order[1] && !is.null(phi.start)) stop("phi.start should have ", 
                order[1], " elements \n")
            phi <- if (is.null(phi.start)) runif(order[1]) else phi.start
            params <- c(params, phi = phi)
            lowerBounds <- c(lowerBounds, phi = rep(-1, order[1]))
            upperBounds <- c(upperBounds, phi = rep(1, order[1]))
        }, {
            if (length(theta.start) != order[2] && !is.null(theta.start)) stop("theta.start should have ", 
                order[2], " elements \n")
            theta <- if (is.null(theta.start)) runif(order[2]) else theta.start
            params <- c(params, theta = theta)
            lowerBounds <- c(lowerBounds, theta = rep(-1, order[2]))
            upperBounds <- c(upperBounds, theta = rep(1, order[2]))
        }, {
            if (length(phi.start) != order[1] && !is.null(phi.start)) stop("phi.start should have ", 
                order[1], " elements \n")
            if (length(theta.start) != order[2] && !is.null(theta.start)) stop("theta.start should have ", 
                order[2], " elements \n")
            phi <- if (is.null(phi.start)) runif(order[1]) else phi.start
            theta <- if (is.null(theta.start)) runif(order[2]) else theta.start
            params <- c(params, phi = phi, theta = theta)
            lowerBounds <- c(lowerBounds, phi = rep(-1, order[1]), 
                theta = rep(-1, order[2]))
            upperBounds <- c(upperBounds, phi = rep(1, order[1]), 
                theta = rep(1, order[2]))
        })
    }
    if ("sigma" %in% names(fam$parameters)) {
        sigma <- fitted(m0, "sigma")[1]
        names(sigma) <- ""
        params <- c(params, sigma = sigma)
        lowerBounds <- c(lowerBounds, sigma = 0.001)
        upperBounds <- c(upperBounds, sigma = 1e+10)
    }
    if ("nu" %in% names(fam$parameters)) {
        nu <- fitted(m0, "nu")[1]
        names(nu) <- ""
        params <- c(params, nu = nu)
        lowerBounds <- c(lowerBounds, nu = -Inf)
        upperBounds <- c(upperBounds, nu = Inf)
    }
    if ("tau" %in% names(fam$parameters)) {
        tau <- fitted(m0, "tau")[1]
        names(tau) <- ""
        params <- c(params, tau = tau)
        lowerBounds <- c(lowerBounds, tau = -Inf)
        upperBounds <- c(upperBounds, tau = Inf)
    }
    fit <- try(nlminb(start = params, objective = garmaLL, lower = lowerBounds, 
        upper = upperBounds, control = control))
    if (any(class(fit) %in% "try-error")) {
        cat("OOPS IT FAILED: let us try once more \n")
        fit <- do.call("garmaFit2", args = as.list(garmacall[-1]))
        return(fit)
    }
    fit1 <- garmaLL(fit$par, save = TRUE)
    df.fit <- length(fit$par)
    Hessian <- hessian(fit$par)
    cat("deviance of  garma model= ", 2 * fit$objective, "\n")
    if (2 * fit$objective > deviance(m0)) 
        warning("There is a problem with the fitted GARMA model here \n")
    out <- list(family = fam$family, parameters = as.character(names(fam$par)), 
        type = fam$type, call = garmacall, y = y, weights = w, 
        G.deviance = 2 * fit$objective, N = N, df.fit = length(fit$par), 
        df.residual = N - df.fit, aic = 2 * fit$objective + 2 * 
            df.fit, sbc = 2 * fit$objective + log(N) * df.fit, 
        method = "nlminb", vcov = solve(Hessian), coef = fit$par)
    if ("mu" %in% names(fam$parameters)) {
        mu <- as.vector(fit1$mu)
        out$mu <- mu
        if (fam$nopar > 1) {
            final <- which(names(fit$par) == "sigma") - 1
            out$mu.coefficients <- fit$par[1:final]
        }
        else {
            out$mu.coefficients <- fit$par
        }
    }
    if ("sigma" %in% names(fam$parameters)) {
        sigma <- fit$par["sigma"]
        out$sigma <- sigma
        out$sigma.coefficients <- fam$sigma.linkinv(sigma)
    }
    if ("nu" %in% names(fam$parameters)) {
        nu <- fit$par["nu"]
        out$nu <- nu
        out$nu.coefficients <- fam$nu.linkinv(nu)
    }
    if ("tau" %in% names(fam$parameters)) {
        tau <- fit$par["tau"]
        out$tau <- tau
        out$tau.coefficients <- fam$tau.linkinv(tau)
    }
    out$residuals <- eval(fam$rqres)
    out$rqres <- fam$rqres
    class(out) <- c("garma", "gamlss")
    out
	}
