\name{bugs2r.distributions}
\alias{bugs2r.distributions}
\title{convert BUGS parameterizations to R paramaterizations}
\usage{
  bugs2r.distributions(..., direction = "bugs2r")
}
\arguments{
  \item{priors}{data.frame with columns distn =
  distribution name, parama, paramb using BUGS default
  parameterizations}
}
\value{
  priors dataframe using R default parameterizations
}
\description{
  convert BUGS parameterizations to R paramaterizations
}
\details{
  bugs2r.distributions is the inverse of
  r2bugs.distributions
}
\examples{
priors <- data.frame(distn = c('weibull', 'lnorm', 'norm', 'gamma'),
                     parama = c(1, 1, 1, 1),
                     paramb = c(2, 2, 2, 2))
bugs2r.distributions(priors)
}
\author{
  David LeBauer, Ben Bolker
}

