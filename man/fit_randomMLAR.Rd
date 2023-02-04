\name{fit_randomMLAR}
\alias{fit_randomMLAR}
\title{Fit the Multilevel AR Model with Random Residual Variances}
\usage{
fit_randomMLAR()
}
\description{
Bayesian estimation of the multilevel autoregressive model. Random effects can be
specified for all model parameters. Note, however, that random residual variances
follow a log-normal distribution.
Additionally, on Level-2, person parameters can be used either as predictors or
as an outcome. That is, an external variable can either be used to predict the
intercepts, the autocorrelation and/or the residual variances, or the
external variable can be predicted by these parameters.
Furthermore, the function optionally returns the model implied mean successive
squared differences of a time series (see: )
}
\examples{
data <- simulate_randomMLAR()
fit_randomMLAR(data)
}
