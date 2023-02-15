\name{fit_mlMSAR}
\alias{fit_mlMSAR}
\title{Fit the Multilevel Markov switching AR Model}
\usage{
fit_mlMSAR()
}
\description{
Bayesian estimation of the multilevel Markov switching AR model. Random effects can be specified for all model parameters. See:
}
\examples{
data <- simulate_mlMSAR()
fit_mlMSAR(data)
}
