##################################################################
## mean and median bias reduction in beta binomial regression ####
##################################################################

## This file contains functions that have been copied
## from package 'betareg' version 3.1-0.
## The authors of betareg are:
## Achim Zeileis <Achim.Zeileis@R-project.org>
## Francisco Cribari-Neto <cribari@ufpe.br>
## Ioannis Kosmidis <ioannis@stats.ucl.ac.uk>
## Alexandre B. Simas 
## Andrea V. Rocha 

##  'brbetabinomial'  were written using as basis the code 
## of 'betareg.

# the functions print.brbetabinomial, summary.brbetabinomial, print.summary.brbetabinomial
# are the modified and adapted  code of  print.betareg, summary.betareg, 
#print.summary.betareg, respectively.    


## Euloge Clovis Kenne Pagui <kenne@stat.unipd.it> [27/07/2019]


#library(VGAM)
#library(Formula)

brbetabinomial <- function(formula,data,weights = NULL, offset = NULL,
                       link.mu="logit",link.phi="identity",start=NULL,type=c("AS_median","AS_mean","AS_ml"),
                       maxit=50,slowit=1,epsilon=1e-08,max_step_factor=11,
                       model = TRUE,y = TRUE, x = FALSE)
{
  ## call
  cl <- match.call()
  if(missing(data)) data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action", "weights", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  
  ## formula
  oformula <- as.formula(formula)
  formula <- as.Formula(formula)
  if(length(formula)[2L] < 2L) {
    formula <- as.Formula(formula(formula), ~ 1)
    simple_formula <- TRUE
  } else {
    if(length(formula)[2L] > 2L) {
      formula <- Formula(formula(formula, rhs = 1:2))
      warning("formula must not have more than two RHS parts")
    }
    simple_formula <- FALSE
  }
  mf$formula <- formula
  
  ## evaluate model.frame
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  
  ## extract terms, model matrix, response
  mt <- terms(formula, data = data)
  mtX <- terms(formula, data = data, rhs = 1L)
  mtZ <- delete.response(terms(formula, data = data, rhs = 2L))
  Y <- model.response(mf, "numeric")
  X <- model.matrix(mtX, mf)
  Z <- model.matrix(mtZ, mf)
  
  # ## sanity checks on response variable (Y)
  # if(length(Y) < 1) stop("empty model")
  # if(!(min(Y) >= 0 & max(Y) <= 1)) stop("invalid dependent variable, all observations must be in (0, 1)")
  # 
  ## convenience variables
  n <- length(Y)
  
  ## type of estimator
  type <- match.arg(type)
  
  ## check links
  if(is.character(link.mu)) link.mu <- match.arg(link.mu,c("logit", "probit", "cloglog", "identity"))
  if(is.character(link.phi)) link.phi <- match.arg(link.phi,c("logit", "probit", "cloglog", "identity"))
  
  ## weights
  weights <- model.weights(mf)
  if(is.null(weights)) weights <- 1
  if(length(weights) == 1) weights <- rep.int(weights, n)
  weights <- as.vector(weights)
  names(weights) <- rownames(mf)
  
  ## offsets
  expand_offset <- function(offset) {
    if(is.null(offset)) offset <- 0
    if(length(offset) == 1) offset <- rep.int(offset, n)
    as.vector(offset)
  }
  ## offsets in mean part of formula
  offsetX <- expand_offset(model.offset(model.part(formula, data = mf, rhs = 1L, terms = TRUE)))
  ## offsets in precision part of formula
  offsetZ <- expand_offset(model.offset(model.part(formula, data = mf, rhs = 2L, terms = TRUE)))
  ## in offset argument (used for mean)
  if(!is.null(cl$offset)) offsetX <- offsetX + expand_offset(mf[, "(offset)"])
  ## collect
  offset <- list(mean = offsetX, precision = offsetZ)
  # print(weights)
  # print(offset)
  # print(Y)
  # print(X)
  # print(Z)
  ## call the actual workhorse: betareg.fit()
  fitmodel <- brbetabinomial.fit( X, Z,Y, weights, offset, link.mu, link.phi,start,
                                 type, maxit,slowit,epsilon,max_step_factor)
  
  ## further model information
  fitmodel$call <- cl
  fitmodel$formula <- oformula
  fitmodel$terms <- list(mean = mtX, precision = mtZ, full = mt)
  fitmodel$levels <- list(mean = .getXlevels(mtX, mf), precision = .getXlevels(mtZ, mf), full = .getXlevels(mt, mf))
  fitmodel$contrasts <- list(mean = attr(X, "contrasts"), precision = attr(Z, "contrasts"))
  if(model) fitmodel$model <- mf
  if(y) fitmodel$y <- Y
  if(x) fitmodel$x <- list(mean = X, precision = Z)
  
  class(fitmodel) <- c("brbetabinomial")
  return(fitmodel)
}

brbetabinomial.fit <- function(x, z, y, weights = NULL, offset = NULL,
                     link.mu="logit",link.phi="identity",start=NULL,type="AS_mean",
                     maxit=50,slowit=1,epsilon=1e-08,max_step_factor=11)
{
  x <- as.matrix(x)
  z <- as.matrix(z)
  n <- NROW(x)
  p <- NCOL(x)
  q <- NCOL(z)
  
  if(is.null(offset)) offset<-list(rep(0,n),rep(0,n))
 
  mu_linkobj <- enrichwith::enrich(make.link(link.mu),with="d2mu.deta")
  mu_linkfun <- mu_linkobj$linkfun
  mu_linkinv <- mu_linkobj$linkinv
  mu_dmu.deta <- mu_linkobj$mu.eta
  mu_d2mu.deta <- mu_linkobj$d2mu.deta
  
  phi_linkobj <- enrichwith::enrich(make.link(link.phi),with="d2mu.deta")
  phi_linkfun <- phi_linkobj$linkfun
  phi_linkinv <- phi_linkobj$linkinv
  phi_dmu.deta <- phi_linkobj$mu.eta
  phi_d2mu.deta <- phi_linkobj$d2mu.deta
  
  ## useful quantities
  Ej <- function(mu,phi,j) {((1-phi)*mu+j*phi)}
  Fj <- function(mu,phi,j) {((1-phi)*(1-mu)+j*phi)}
  Gj <- function(phi,j) {((1-phi)+j*phi)}
  
  ### required quantities ##
  
  key_quantities <- function(pars)
  {
    betas <- pars[1:p]
    gammas <- pars[-c(1:p)]
    mu_etas <- drop(x %*% betas + offset[[1L]])
    mus <- mu_linkinv(mu_etas)
    d1mus <- mu_dmu.deta(mu_etas)
    d2mus <- mu_d2mu.deta(mu_etas)
    phi_etas <- drop(z %*% gammas + offset[[2L]])
    phis <- phi_linkinv(phi_etas)
    d1phis <- phi_dmu.deta(phi_etas)
    d2phis <- phi_d2mu.deta(phi_etas)
    # dens<-lapply(seq_along(weights), function(t) dbetabinom(x=0:weights[t], 
    #       size=weights[t], prob=mus[t], rho=phis[t])) ## density evaluated in 0:m_i
    E1<-E2<-E3<-E4<-E5<-E6<-E7<-E8<-rep(NA,n)
    E9<-E10<-E11<-E12<-E13<-E14<-E15<-E16<-rep(NA,n)
    for(i in  seq_along(weights))
    {
      index <- 0:weights[i]
      dens <- dbetabinom(x=index,size=weights[i], prob=mus[i], rho=phis[i])
      l_mu <- scoremu(mus[i],phis[i],index,weights[i])
      l_phi <- scorephi(mus[i],phis[i],index,weights[i])
      l_mumu <- scoremumu(mus[i],phis[i],index,weights[i])
      l_muphi <- scoremuphi(mus[i],phis[i],index,weights[i])
      l_phiphi <- scorephiphi(mus[i],phis[i],index,weights[i])
      E1[i] <- sum(l_mumu*dens)
      E2[i] <- sum(l_muphi*dens)
      E3[i] <- sum(l_phiphi*dens)
      E4[i] <- sum(l_mu^3*dens)
      E5[i] <- sum(l_mu^2*l_phi*dens)
      E6[i] <- sum(l_mu*l_phi^2*dens)
      E7[i] <- sum(l_phi^3*dens)
      E8[i] <- sum(l_mumu*l_mu*dens)
      E9[i] <- sum(l_mu^2*dens)
      E10[i] <- sum(l_mu*l_muphi*dens)
      E11[i] <- sum(l_mu*l_phiphi*dens)
      E12[i] <- sum(l_mu*l_phi*dens)
      E13[i] <- sum(l_phi*l_mumu*dens)
      E14[i] <- sum(l_phi*l_muphi*dens)
      E15[i] <- sum(l_phi*l_phiphi*dens)
      E16[i] <- sum(l_phi^2*dens)
    }
    out <- list(betas=betas, gammas=gammas, mu_etas=mu_etas, 
                phi_etas=phi_etas, mus=mus, phis=phis,d1mus=d1mus,
                d2mus=d2mus, d1phis=d1phis, d2phis=d2phis,E1=E1,E2=E2,E3=E3,
                E4=E4,E5=E5,E6=E6,E7=E7,E8=E8,E9=E9,E10=E10,E11=E11,E12=E12,
                E13=E13,E14=E14,E15=E15,E16=E16)
  }
  
  ## score function for mu and phi (only one observation)##
  
  scoremu <- Vectorize( function(mu, phi,y,m){ ## m=weights
    sum1 <- sum2 <- 0
    
    if(y!=0)
    {
      index <- 0:(y-1)
      sum1 <-  sum(1/Ej(mu,phi,index))
    }
    
    if(y<m)
    {
      index <- 0:(m-y-1)
      sum2 <-  sum(1/Fj(mu,phi,index))
    }
    
    (1-phi)*(sum1-sum2)
    
  }, "y")
  
  
  scorephi <- Vectorize( function(mu, phi,y,m){ ## m=weights
    sum1 <- sum2 <- sum3 <- 0
    
    if(y!=0)
    {
      index <- 0:(y-1)
      sum1 <- sum((index-mu)/Ej(mu,phi,index))
    }
    
    if(y<m)
    {
      index <- 0:(m-y-1)
      sum2 <- sum((index+mu-1)/Fj(mu,phi,index))
    }
    
    index <- 0:(m-1)
    sum3 =  sum((index-1)/Gj(phi,index))
    rm(index)
    (sum1 + sum2 - sum3)
  }, "y")
  
  ### second derivatives for one observation  ##

  
  scoremumu <- Vectorize( function(mu, phi, y, m){ ## m=weights
    
    sum1 <- sum2 <- 0
    
    if(y!=0)
    {
      index <- 0:(y-1)
      sum1 <-  sum(1/Ej(mu,phi,index)^2)
    }
    
    if(y<m)
    {
      index <- 0:(m-y-1)
      sum2 <-  sum(1/Fj(mu,phi,index)^2)
    }
    
    -(1-phi)^2*(sum1+sum2)
    
  }, "y")
  
  scoremuphi <- Vectorize( function(mu, phi, y, m){ ## m=weights
    
    sum1 <- sum2 <- 0
    
    if(y!=0)
    {
      index <- 0:(y-1)
      sum1 <-  sum(index/Ej(mu,phi,index)^2)
    }
    
    if(y<m)
    {
      index <- 0:(m-y-1)
      sum2 <-  sum(index/Fj(mu,phi,index)^2)
    }
    
    -(sum1-sum2)
    
  }, "y")
  
  scorephiphi <- Vectorize( function(mu, phi, y, m){ ## m=weights
    sum1 <- sum2 <- sum3 <- 0
    
    if(y!=0)
    {
      index <- 0:(y-1)
      sum1 <- sum((mu-index)^2/Ej(mu,phi,index)^2)
    }
    
    if(y<m)
    {
      index <- 0:(m-y-1)
      sum2 <- sum((index+mu-1)^2/Fj(mu,phi,index)^2)
    }
    
    index <- 0:(m-1)
    sum3 =  sum((index-1)^2/Gj(phi,index)^2)
    rm(index)
    -(sum1 + sum2 - sum3)
  }, "y")
  
  ## negative likelihood ##
  nllik <- function(pars,fit=NULL)
  {
    if (is.null(fit))
    {
      fit <- key_quantities(pars)
    }
    with(fit,{
      # print(mus)
      # print(phis)
      # print(weights)
      -sum(dbetabinom(x=y,size=weights, prob=mus, rho=phis,log = TRUE))
    })
  }
  
  ## score function##
  
  score <- function(pars,fit=NULL)
  {
    if (is.null(fit))
    {
      fit <- key_quantities(pars)
    }
    with(fit,{
     colSums(cbind( x*d1mus*sapply(1:n,function(i)scoremu(mus[i],phis[i],y[i],weights[i])),
                z*d1phis*sapply(1:n,function(i)scorephi(mus[i],phis[i],y[i],weights[i]))))
    })
  }
  
  ## Fisher information
  
  information <- function(pars,inverse=FALSE, fit= NULL)
  {
    if (is.null(fit))
    {
      fit <- key_quantities(pars)
    }
    with(fit,{
      wbb <- -E1*d1mus^2 
      ibb <- crossprod(wbb*x,x)
      wbg <- -E2*d1mus*d1phis
      ibg <- crossprod(wbg*x,z)
      wgg <- -E3*d1phis^2
      igg <- crossprod(wgg*z,z)
      info <- rbind(cbind(ibb,ibg),cbind(t(ibg),igg))
      if (!inverse) 
        info
      else chol2inv(chol(info))
    })
  }
  
  ## mean  adjustment term
  
  mean_adjustment <- function(pars,fit= NULL)
  {
    if (is.null(fit))
    {
      fit <- key_quantities(pars)
    }
    invInfo <- information(pars,TRUE,fit)
    with(fit,{
     PQsum <- function(t)
     {
       if ( t <= p )
       {
         Xt <- x[,t]
         bb <- crossprod(x,((E4+E8)*d1mus^3+E9*d1mus*d2mus)*Xt*x)
         bg <- crossprod(x,(E5+E10)*d1mus^2*d1phis*Xt*z)
         gg <- crossprod(z,((E6+E11)*d1mus*d1phis^2+E12*d1mus*d2phis)*Xt*z)
       }
       else
       {
         Zt <- z[, t-p]
         bb <- crossprod(x,((E5+E13)*d1mus^2*d1phis+E12*d1phis*d2mus)*Zt*x)
         bg <- crossprod(x,(E6+E14)*d1mus*d1phis^2*Zt*z)
         gg <- crossprod(z,((E7+E15)*d1phis^3+E16*d1phis*d2phis)*Zt*z)
       }
       pq<-rbind(cbind(bb,bg),cbind(t(bg),gg))
       sum(diag(invInfo %*% pq))/2
     }
     sapply(1:(p + q), PQsum)
    })
  } 
  
  
  ## median  adjustment term
  
  median_adjustment <- function(pars,fit= NULL)
  {
    if (is.null(fit))
    {
      fit <- key_quantities(pars)
    }
    invInfo <- information(pars,TRUE,fit)
    info <- information(pars,FALSE,fit)
    m1_vector <- numeric((p+q))
    with(fit,{
      PQsum <- function(t)
      {
        if ( t <= p )
        {
          Xt <- x[,t]
          bb <- crossprod(x,((E4+E8)*d1mus^3+E9*d1mus*d2mus)*Xt*x)
          bg <- crossprod(x,(E5+E10)*d1mus^2*d1phis*Xt*z)
          gg <- crossprod(z,((E6+E11)*d1mus*d1phis^2+E12*d1mus*d2phis)*Xt*z)
        }
        else
        {
          Zt <- z[, t-p]
          bb <- crossprod(x,((E5+E13)*d1mus^2*d1phis+E12*d1phis*d2mus)*Zt*x)
          bg <- crossprod(x,(E6+E14)*d1mus*d1phis^2*Zt*z)
          gg <- crossprod(z,((E7+E15)*d1phis^3+E16*d1phis*d2phis)*Zt*z)
        }
        pq<-rbind(cbind(bb,bg),cbind(t(bg),gg))
        sum(diag(invInfo %*% pq))/2
      }
      
      PQsum2 <- function(t)
      {
        if ( t <= p )
        {
          Xt <- x[,t]
          bb <- crossprod(x,((E4/3+E8/2)*d1mus^3+E9*d1mus*d2mus/2)*Xt*x)
          bg <- crossprod(x,(E5/3+E10/2)*d1mus^2*d1phis*Xt*z)
          gg <- crossprod(z,((E6/3+E11/2)*d1mus*d1phis^2+E12*d1mus*d2phis/2)*Xt*z)
        }
        else
        {
          Zt <- z[, t-p]
          bb <- crossprod(x,((E5/3+E13/2)*d1mus^2*d1phis+E12*d1phis*d2mus/2)*Zt*x)
          bg <- crossprod(x,(E6/3+E14/2)*d1mus*d1phis^2*Zt*z)
          gg <- crossprod(z,((E7/3+E15/2)*d1phis^3+E16*d1phis*d2phis/2)*Zt*z)
        }
        pq<-rbind(cbind(bb,bg),cbind(t(bg),gg))
        sum(diag(h_r %*% pq))
      }
      for(r in 1:(p + q))
      {
        InvInfo_r <- invInfo[r, ]
        h_r <- tcrossprod(InvInfo_r)/InvInfo_r[r]
        m1_vector[r] <- sum(InvInfo_r*sapply(1:(p + q),PQsum2))
      }
      sapply(1:(p + q), PQsum)-drop(info %*% m1_vector) 
    })
  } 
  
  ## adjustment term function ##
  
  adjustment_function <- switch(type, 
                                AS_mean = mean_adjustment, 
                                AS_median = median_adjustment,
                                AS_ml = function(pars, ...) 0)
  
  ## required components for computing the adjusted scores ##
  
  compute_step_components <- function(pars,  fit = NULL)
  {
    if (is.null(fit))
    {
      fit <- key_quantities(pars)
    }
    grad <- score(pars,fit)
    invInfo <- try(information(pars, inverse = TRUE,fit = fit))
    info <- try(information(pars, inverse = FALSE,fit = fit))
    failed_inversion <- inherits(invInfo, "try-error")
    adjustment <- try(adjustment_function(pars,fit))
    failed_adjustment <- inherits(adjustment, "try-error")
    out <- list(grad = grad, info=info, invInfo = invInfo, 
        adjustment = adjustment, failed_inversion = failed_inversion,
        failed_adjustment=failed_adjustment)
    out
  }  
  
  adjusted_grad <- function(pars)
  {
    fit <- key_quantities(pars)
    components <- compute_step_components(pars,fit)
    with(components, 
         {
           grad + adjustment
         })
  }
 # matse <- NULL 
  par <- start
  quantities <- key_quantities(par)
  step_components <- compute_step_components(start, fit = quantities)
  if (step_components$failed_inversion) {
    warning("failed to invert the information matrix")
  }
  if (step_components$failed_adjustment) {
    warning("failed to calculate score adjustment")
  }
  adjusted_grad <- with(step_components, 
                        {
                          grad + adjustment
                        })
  step_par <- drop(step_components$invInfo %*% adjusted_grad)
  
  for (iter in seq.int(maxit)) 
  {
    step_factor <- 0
    testhalf <- TRUE
    while (testhalf & step_factor < max_step_factor) 
    {
      step_par_previous <- step_par
      par <- par + slowit*2^(-step_factor)*step_par
      
      quantities <- key_quantities(par)
      step_components <- compute_step_components(start, fit = quantities)
      if (step_components$failed_inversion) 
      {
        warning("failed to invert the information matrix")
      }
      if (step_components$failed_adjustment) {
        warning("failed to calculate score adjustment")
      }
      adjusted_grad <- with(step_components, 
                            {
                              grad + adjustment
                            })
      step_par <- drop(step_components$invInfo %*% adjusted_grad)
      
      if (step_factor == 0 & iter == 1) 
      {
        testhalf <- TRUE
      }
      else 
      {
        testhalf <- sum(abs(step_par),na.rm=TRUE) > sum(abs(step_par_previous),na.rm = TRUE)
      }
      step_factor <- step_factor + 1
      # if (control$trace) {
      #   trace_iteration()
      # }
    }
    #if(type=="AS_mean")
    #{
      #matse <- rbind(matse,sqrt(diag(step_components$invInfo)))
    #}
    failed <- step_components$failed_inversion | step_components$failed_adjustment
    if (failed | sum(abs(step_par),na.rm = TRUE) < epsilon) {
      break
    }
  }
  
  
  if(iter >= maxit)
  {
    convergence <- 0
    warning("optimization failed to converge")
  }
  else
  {
    convergence <- 1
  }
  
  # if(type=="AS_ML")
  # {
  #   if(link.mu=="identity" & link.phi=="identity")
  #   ml<- nlminb(start,nllik,lower = rep(1e-04,p+q),upper = rep(0.99,p+q))
  #   if(link.mu!="identity" & link.phi=="identity")
  #     ml<- nlminb(start,nllik,lower = c(rep(-Inf,p),rep(1e-04,q)),upper = c(rep(+Inf,p),rep(0.99,q)))
  #   if(link.mu=="identity" & link.phi!="identity")
  #     ml<- nlminb(start,nllik,lower = c(rep(1e-04,p),rep(-Inf,q)),upper = c(rep(0.99,p),rep(+Inf,q)))
  #   if(link.mu!="identity" & link.phi!="identity")
  #     ml<- nlminb(start,nllik)
  #   par<-ml$par 
  #   se <- sqrt(diag(information(mle,TRUE,NULL)))
  #   convergence <- 1-ml$convergence
  #   iter <- ml$iterations
  # }
  
  
  #fit <- key_quantities(par)
  beta <- quantities$betas
  gamma <- quantities$gammas
  vcov <- step_components$invInfo
  #se <- sqrt(diag( vcov))
  ## coefficient names
  names(beta) <-  if(link.mu == "identity") "(mu)" else colnames(x)
  names(gamma) <- if(link.phi == "identity") "(phi)" else colnames(z)
  rownames(vcov) <- colnames(vcov)  <- c(if(link.mu == "identity") "(mu)" else colnames(x),
                                         if(link.phi == "identity") "(phi)" else paste("(phi)", colnames(z), sep = "_"))
  
  
  
  list(coefficients=list(mean = beta, precision = gamma),convergence=convergence,
       iter=iter,typeestimator=type,vcov=vcov,grad=adjusted_grad,
       link = list(mean = link.mu, precision = link.phi))
  #list(coefficients=list(mean = beta, precision = gamma),se=se,conv=convergence,iter=iter,typeEstimator=type)
  #list(coef=par,se=se,conv=convergence,iter=iter)
}  


print.brbetabinomial <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")

  if(!x$convergence) {
    cat("model did not converge\n")
  } else {
    if(length(x$coefficients$mean)) {
      cat(paste("Coefficients (mean model with ", x$link$mean, " link):\n", sep = ""))
      print.default(format(x$coefficients$mean, digits = digits), print.gap = 2, quote = FALSE)
      cat("\n")
    } else cat("No coefficients (in mean model)\n\n")
    
    if(length(x$coefficients$precision)) {
      cat(paste("Phi coefficients (precision model with ", x$link$precision, " link):\n", sep = ""))
      print.default(format(x$coefficients$precision, digits = digits), print.gap = 2, quote = FALSE)
      cat("\n")
    } else cat("No coefficients (in precision model)\n\n")
  
  }
}

summary.brbetabinomial <- function(object,   ...)
{

  ## coefficient table
  p <- length(object$coefficients$mean)
  q <- length(object$coefficients$precision)
  cf <- as.vector(do.call("c", object$coefficients))
  se <- sqrt(diag(object$vcov))
  cf <- cbind(cf, se, cf/se, 2 * pnorm(-abs(cf/se)))
  colnames(cf) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  cf <- list(mean = cf[seq.int(length.out = p), , drop = FALSE], precision = cf[seq.int(length.out = q) + p, , drop = FALSE])
  rownames(cf$mean) <- names(object$coefficients$mean)
  rownames(cf$precision) <- names(object$coefficients$precision)
  object$coefficients <- cf

  ## number of iterations
  
  object$iterations <- object$iter

  ## delete some slots
  object$terms <- object$model <- object$y <-NULL
    object$x <- object$levels <- object$contrasts <-  NULL

  ## return
  class(object) <- "summary.brbetabinomial"
  object
}

print.summary.brbetabinomial <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")
  
  if(!x$convergence) {
    cat("model did not converge\n")
  }
    
  if(NROW(x$coefficients$mean)) {
    cat(paste("\nCoefficients (mean model with ", x$link$mean, " link):\n", sep = ""))
    printCoefmat(x$coefficients$mean, digits = digits, signif.legend = FALSE)
  } else cat("\nNo coefficients (in mean model)\n")
    
    
  if(NROW(x$coefficients$precision)) {
    cat(paste("\nPhi coefficients (precision model with ", x$link$precision, " link):\n", sep = ""))
    printCoefmat(x$coefficients$precision, digits = digits, signif.legend = FALSE)
  } else cat("\nNo coefficients (in precision model)\n")
    
    
    if(getOption("show.signif.stars") & any(do.call("rbind", x$coefficients)[, 4L] < 0.1))
      cat("---\nSignif. codes: ", "0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1", "\n")
    
    cat("\nType of estimator:", x$typeestimator, switch(x$typeestimator,
                                               "AS_ml" = "(maximum likelihood)",
                                               "AS_mean" = "(mean bias-reduced)",
                                               "AS_median" = "(median bias-reduced)"))
    # cat("\nLog-likelihood:", formatC(x$loglik, digits = digits),
    #     "on", sum(sapply(x$coefficients, NROW)), "Df")
    
  cat(paste("\nNumber of iterations in the quasi-Fisher scoring:", x$iterations, "\n"))
    
  
  invisible(x)
}
