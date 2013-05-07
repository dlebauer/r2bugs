\name{bugs2r.distributions}
\alias{bugs2r.distributions}
\title{convert BUGS parameterizations to R paramaterizations}
\usage{
  bugs2r.distributions(..., direction = "bugs2r")
}
\arguments{
  \item{...}{values passed to
  \code{\link{r2bugs.distributions}}, currently a
  data.frame called \code{priors} with columns distn =
  distribution name, parama, paramb using BUGS default
  parameterizations (see example)}

  \item{direction,}{one of 'bugs2r' (default) or 'r2bugs'.
  Determines if function should convert from BUGS to R, or
  vice-versa. The \code{\link{r2bugs.distributions}} does
  the reverse.}
}
\value{
  priors dataframe using R default parameterizations
}
\description{
  convert BUGS parameterizations to R paramaterizations
}
\details{
  \code{\link{bugs2r.distributions}} is the inverse of the
  \code{\link{r2bugs.distributions}}. The only difference
  in these functions is the defalut value of the
  \code{direction} parameter.
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

