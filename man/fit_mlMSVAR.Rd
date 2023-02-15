\name{fit_mlMSVAR}
\alias{fit_mlMSVAR}
\title{Fit the Multilevel Markov switching VAR Model}
\usage{
fit_mlMSVAR()
}
\description{
Bayesian estimation of the multilevel Markov switching vector autoregressive model. Random effects can be
specified for all model parameters.
}
\examples{
data <- simulate_mlMSVAR()
fit_mlMSVAR(data)
}
