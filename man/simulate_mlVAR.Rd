\name{simulate_mlVAR}
\alias{simulate_mlVAR}
\title{Simulate data under the Multilevel VAR Model}
\usage{
simulate_mlVAR()
}
\description{
Simulate data under the multilevel vector autoregressive model. Random effects can be
simulated for all model parameters.
}
\examples{
data <- simulate_mlVAR()
fit_mlVAR(data)
}
