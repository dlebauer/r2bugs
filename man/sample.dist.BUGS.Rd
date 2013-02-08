\name{sample.dist.BUGS}
\alias{sample.dist.BUGS}
\title{Sample from an R distribution using JAGS}
\usage{
  sample.dist.BUGS(prior = data.frame(distn = "norm", parama = 0, paramb = 1),
    n.iter = 1e+05, n = NULL)
}
\arguments{
  \item{prior}{dataframe with distribution name and
  parameters}

  \item{n.iter}{number of mcmc samples}

  \item{n}{number of samples returned, if NULL (default),
  output will have n.iter/4 samples}
}
\value{
  vector of samples
}
\description{
  Use JAGS to sample from a distribution stated using R
  parameterization
}
\details{
  Takes a distribution with R parameterization, converts it
  to a BUGS parameterization, and then samples from the
  distribution using JAGS
}
\author{
  David LeBauer
}

