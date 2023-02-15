\name{fit_mlVAR}
\alias{fit_mlVAR}
\title{Fit the Multilevel VAR Model}
\usage{
fit_mlVAR()
}
\description{
Bayesian estimation of the multilevel vector autoregressive model. Random effects can be
specified for all model parameters.
}
\examples{
data <- simulate_mlVAR()
fit_mlVAR(data)
}
