#-------------------------------------------------------------------------------
# Copyright (c) 2012 University of Illinois, NCSA.
# All rights reserved. This program and the accompanying materials
# are made available under the terms of the 
# University of Illinois/NCSA Open Source License
# which accompanies this distribution, and is available at
# http://opensource.ncsa.illinois.edu/license.html
#-------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------#
##' convert R parameterizations to BUGS paramaterizations
##' 
##' R and BUGS have different parameterizations for some distributions. This function transforms the distributions from R defaults to BUGS defaults. BUGS is an implementation of the BUGS language, and these transformations are expected to work for bugs.
##' @title convert R parameterizations to BUGS paramaterizations
##' @param priors data.frame with columns distn = distribution name, parama, paramb using R default parameterizations
##' @param direction, one of 'r2bugs' (default), 'bugs2r'. Determines if function should convert from R to BUGS, or vice-versa. The \code{\link{bugs2r.distributions}} does the reverse.
##' @return priors dataframe using JAGS default parameterizations
##' @author David LeBauer, Ben Bolker
##' @export
##' @examples
##' priors <- data.frame(distn = c('weibull', 'lnorm', 'norm', 'gamma'),
##'                      parama = c(1, 1, 1, 1),
##'                      paramb = c(2, 2, 2, 2))
##' r2bugs.distributions(priors)
##' 

r2bugs.distributions <- function(priors, direction = 'r2bugs') {

  priors$distn  <- as.character(priors$distn)
  priors$parama <- as.numeric(priors$parama)
  priors$paramb <- as.numeric(priors$paramb)

  ## index dataframe according to distribution
  norm   <- priors$distn %in% c('norm', 'lnorm')    # these have same tramsform
  weib   <- grepl("weib", priors$distn)             # matches r and bugs version
  gamma  <- priors$distn == 'gamma'
  chsq   <- grepl("chisq", priors$distn)            # matches r and bugs version
  bin    <- priors$distn %in% c('binom', 'bin')     # matches r and bugs version
  nbin   <- priors$distn %in% c('nbinom', 'negbin') # matches r and bugs version

  ## Check that no rows are categorized into two distributions
  if(max(rowSums(cbind(norm, weib, gamma, chsq, bin, nbin))) > 1) {
    badrow <- rowSums(cbind(norm, weib, gamma, chsq, bin, nbin)) > 1
    stop(paste(unique(priors$distn[badrow])),
         "are identified as > 1 distribution")
  }

  exponent <- ifelse(direction == "r2bugs", -2, -0.5) 
  ## Convert sd to precision for norm & lnorm
  priors$paramb[norm] <-  priors$paramb[norm] ^ exponent
  if(direction == 'r2bugs'){
    ## Convert R parameter b to BUGS parameter lambda by l = (1/b)^a
    priors$paramb[weib] <-   (1 / priors$paramb[weib]) ^ priors$parama[weib]
  } else if (direction == 'bugs2r') {
    ## Convert BUGS parameter lambda to BUGS parameter b by b = l^(-1/a)
    priors$paramb[weib] <-  priors$paramb[weib] ^ (- 1 / priors$parama[weib] )
 
  }
  ## Reverse parameter order for binomial and negative binomial
  priors[bin | nbin, c('parama', 'paramb')] <-  priors[bin | nbin, c('paramb', 'parama')]
  
  ## Translate distribution names
  if(direction == "r2bugs"){
    priors$distn[weib] <- "weib"
    priors$distn[chsq] <- "chisqr"
    priors$distn[bin]  <- "bin"
    priors$distn[nbin] <- "negbin"
  } else if(direction == "bugs2r"){
    priors$distn[weib] <- "weibull"
    priors$distn[chsq] <- "chisq"
    priors$distn[bin]  <- "binom"
    priors$distn[nbin] <- "nbinom"
  }
  return(priors)
}

##' convert BUGS parameterizations to R paramaterizations
##' 
##' \code{\link{bugs2r.distributions}} is the inverse of the \code{\link{r2bugs.distributions}}. The only difference in these functions is the defalut value of the \code{direction} parameter.
##' @title convert BUGS parameterizations to R paramaterizations
##' @param ... values passed to \code{\link{r2bugs.distributions}}, currently a data.frame called \code{priors} with columns distn = distribution name, parama, paramb using BUGS default parameterizations (see example)
##' @param direction, one of 'bugs2r' (default) or 'r2bugs'. Determines if function should convert from BUGS to R, or vice-versa. The \code{\link{r2bugs.distributions}} does the reverse.
##' @return priors dataframe using R default parameterizations
##' @author David LeBauer, Ben Bolker
##' @export
##' @examples
##' priors <- data.frame(distn = c('weibull', 'lnorm', 'norm', 'gamma'),
##'                      parama = c(1, 1, 1, 1),
##'                      paramb = c(2, 2, 2, 2))
##' bugs2r.distributions(priors)
bugs2r.distributions <- function(..., direction = "bugs2r") {
  return(r2bugs.distributions(..., direction))
}

##' Use JAGS to sample from a distribution stated using R parameterization
##'
##' Takes a distribution with R parameterization, converts it to a
##' BUGS parameterization, and then samples from the distribution using
##' JAGS
##' @title Sample from an R distribution using JAGS
##' @param prior dataframe with distribution name and parameters 
##' @param n.iter number of mcmc samples
##' @param n number of samples returned, if NULL (default),
##' output will have n.iter/4 samples 
##' @return vector of samples
##' @export
##' @author David LeBauer
sample.dist.BUGS <- function(prior = data.frame(
                         distn = "norm",
                         parama = 0,
                         paramb = 1),
                       n.iter = 100000,
                       n = NULL) {
  require(rjags)
  model.string <- paste0("model{Y ~ d", prior$distn, "(",
                         prior$parama, 
                         ifelse(prior$distn == "chisq", "", paste0(", ", prior$paramb)),
                         ")\n a <- x}")
  
  writeLines(model.string, con = "test.bug")
  j.model  <- jags.model(file = "test.bug", data = list(x = 1))
  mcmc.object <- window(coda.samples(model = j.model,
                              variable.names = c('Y'),
                              n.iter = n.iter*4,
                              thin = 2),
                        start = n.iter / 2)
  Y <- as.matrix(mcmc.object)[,"Y"]
  if(is.null(n)){
    n <- n.iter
  }
  Y <- sample(Y, n)
  return(Y)
}

##' Take n random samples from prior
##'
##' @title Sample from a probability distirbution using R
##' @param prior data.frame with distn, parama, paramb
##' @param n number of samples to return
##' @return vector with n random samples from prior
##' @seealso \link{pr.samp}
##' @export
#--------------------------------------------------------------------------------------------------#
sample.dist.R <- function(prior, n) {
  do.call(paste('r', prior$distn, sep=""), list(n, prior$parama, prior$paramb))
}
