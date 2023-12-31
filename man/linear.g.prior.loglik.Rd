% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/likelihoods.R
\name{linear.g.prior.loglik}
\alias{linear.g.prior.loglik}
\title{Log likelihood function for linear regression using Zellners g-prior}
\usage{
linear.g.prior.loglik(y, x, model, complex, params = list(g = 4))
}
\arguments{
\item{y}{A vector containing the dependent variable}

\item{x}{The matrix containing the precalculated features}

\item{model}{The model to estimate as a logical vector}

\item{complex}{A list of complexity measures for the features}

\item{params}{A list of parameters for the log likelihood, supplied by the user}
}
\value{
A list with the log marginal likelihood combined with the log prior (crit) and the posterior mode of the coefficients (coefs).
}
\description{
Log likelihood function for linear regression using Zellners g-prior
}
\examples{
linear.g.prior.loglik(rnorm(100), matrix(rnorm(100)), TRUE, list(oc=1))

}
