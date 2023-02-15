\name{fit_mlHMM}
\alias{fit_mlHMM}
\title{fit the Multilevel Hidden Markov Model with a multivariate normal conditional distribution}
\usage{
fit_mlHMM()
}
\description{
fit the Multilevel Hidden Markov Model with a multivariate normal conditional distribution. The function allows to specify random effects in all model parameters. See details.
}
\examples{
data <- simulate_mlHMM()
fit_mlHMM(data)
}
