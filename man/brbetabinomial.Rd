\name{brbetabinomial}
\alias{brbetabinomial}
\alias{brbetabinomial.fit}
\alias{print.brbetabinomial}
\alias{summary.brbetabinomial}
\alias{print.summary.brbetabinomial}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Mean and median bias reduction in beta-binomial regression
}
\description{
Fit beta-binomial regression using the adjusted score functions to  median bias reduction (Kenne Pagui et al., 2017) and
mean bias reduction (Firth, 1993). Both approaches are effective in solving the boundary estimate problem that may arise for maximum likelihood in beta-binomial regression. 
}
\usage{
brbetabinomial(formula, data, weights = NULL, offset = NULL, link.mu = "logit", link.phi = "identity", start = NULL, type = c("AS_median", "AS_mean", "AS_ml"), maxit = 50, slowit = 1, epsilon = 1e-08, max_step_factor = 11, model = TRUE, y = TRUE, x = FALSE)

brbetabinomial.fit(x, z, y, weights = NULL, offset = NULL,
                     link.mu="logit",link.phi="identity",start=NULL,type="AS_mean",
                     maxit=50,slowit=1,epsilon=1e-08,max_step_factor=11)
}

\arguments{
  \item{formula}{
symbolic description of the model (of type \code{y ~ x}
    or \code{y ~ x | z}; for details see below).
}
  \item{data}{
an optional data frame, list or environment (or object coercible by \code{\link{as.data.frame}} to a data frame) containing the variables in the model. If not found in data, the variables are taken from \code{environment(formula)}. See also \code{\link{glm}}
}
  \item{weights}{
binomial totals.
}
  \item{offset}{
an optional vector specifying a subset of observations to be used in the fitting process.
}
  \item{link.mu}{
link function for mean model.
}
  \item{link.phi}{
link function for precision model.
}
  \item{start}{
starting values.
}
  \item{type}{
type of estimator.
}
  \item{maxit}{
maximum number of iterations needed for the fisher scoring optimization.
}
  \item{slowit}{
a positive real taking values in (0,1] used as a multiplier for the stepsize. 
}

\item{epsilon}{
tollerance for the convergence.
}

\item{max_step_factor}{
maximum number for the step halving.
}

\item{model}{
logical. If \code{TRUE} the model frame is returned.
}

\item{y}{
for \code{brbetabinomial}: \code{y} is a logical value indicates whether the response vector used in the fitting process should be returned as component 
of the returned value

for \code{brbetabinomial.fit}: \code{y} is a vector of observations, representing the 
number of successes.
}
\item{x}{
for \code{brbetabinomial}: \code{x} is a logical value indicates whether the model matrix used in the fitting process should be returned as component 
of the returned value

for \code{brbetabinomial.fit}: \code{x} is a design matrix for the mean model.
}

\item{z}{
numeric matrix. Regressor matrix for the precision model, the default corresponding to the intercept.
}
}
\details{
to be filled.
}

\value{
to be filled.
}

\references{
Firth D. (1993). Bias reduction of maximum likelihood estimates,
Biometrika, **80**, 27-38

Kenne Pagui, E. C., A. Salvan, and N. Sartori (2019). Efficient implementation of median bias reduction.

Kenne Pagui, E. C., A. Salvan, and N. Sartori (2017). Median bias
reduction of maximum likelihood estimates. *Biometrika*, **104**,
923-938
}

\author{
 Euloge Clovis Kenne Pagui[aut, cre] \email{kenne@stat.unipd.it}, 
 Nicola Sartori[aut] \email{sartori@stat.unipd.it},
 Alessandra Salvan[aut] \email{salvan@stat.unipd.it}
}
\note{
to be filled.
}

\seealso{
to be filled.
}
\examples{
\dontrun{
## The results below provide Table 3 of Kenne Pagui et. al (2019)
require(VGAM)
data("lirat",package="VGAM")
## subset of lirat where the total binomials are less than or equal to 11
liratsub=subset(lirat,N<=11)
liratsub$grp=as.factor(liratsub$grp)
mle.vglm.sub <- suppressWarnings(vglm(cbind(R,N-R)~grp+hb, betabinomial(zero = 2,irho=0.2, lrho = "identitylink"),
                data=liratsub ))
mle.brbetabinomial.sub=brbetabinomial(R~grp+hb,weights=N,type="AS_ml",start=c(rep(0,5),0.12)
             ,link.mu = "logit",link.phi = "identity",data=liratsub)                
start.sub = as.vector(do.call("c", mle.brbetabinomial.sub$coefficients))
meanBR.sub=brbetabinomial(R~grp+hb,weights=N,type="AS_mean",start=start.sub
             ,link.mu = "logit",link.phi = "identity",data=liratsub) 
medianBR.sub=brbetabinomial(R~grp+hb,weights=N,type="AS_median",start=start.sub
             ,link.mu = "logit",link.phi = "identity",data=liratsub) 

## whole dataset lirat
lirat$grp <- factor(lirat$grp)

mle.vglm <- suppressWarnings(vglm(cbind(R,N-R)~grp+hb, betabinomial(zero = 2,irho=0.2, lrho = "identitylink"),
                data=lirat ))
mle.brbetabinomial=brbetabinomial(R~grp+hb,weights=N,type="AS_ml",start=c(rep(0,5),0.12)
             ,link.mu = "logit",link.phi = "identity",data=lirat)                
start = as.vector(do.call("c", mle.brbetabinomial$coefficients))
meanBR=brbetabinomial(R~grp+hb,weights=N,type="AS_mean",start=start
             ,link.mu = "logit",link.phi = "identity",data=lirat) 
medianBR=brbetabinomial(R~grp+hb,weights=N,type="AS_median",start=start
             ,link.mu = "logit",link.phi = "identity",data=lirat)  
}             
}

