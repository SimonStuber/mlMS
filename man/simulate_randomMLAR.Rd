\name{simulate_randomMLAR}
\alias{simulate_randomMLAR}
\title{Simulate data under the Multilevel random AR Model}
\usage{
simulate_randomMLAR()
}
\description{
Simulate data under the Multilevel random AR model with random effects in all mode parameters. Additionally, on Level-2, an external variable can be simulated either as an outcome or an predictor of person specific parameters.
}
\examples{
data <- simulate_randomMLAR()
fit_randomMLAR(data)
}
