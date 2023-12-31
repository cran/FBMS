% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/results.R
\name{string.population.models}
\alias{string.population.models}
\title{Function to get a character respresentation of a list of models}
\usage{
string.population.models(features, models, round = 2, link = "I")
}
\arguments{
\item{features}{A list of feature objects on which the models are build}

\item{models}{A list of model objects}

\item{round}{Rounding precision for parameters of the features}

\item{link}{The link function to use, as a string}
}
\value{
A matrix of character representations of a list of models.
}
\description{
Function to get a character respresentation of a list of models
}
\examples{
result <- gmjmcmc(matrix(rnorm(600), 100), P = 2, gaussian.loglik, NULL, c("p0", "exp_dbl"))
string.population.models(result$populations[[2]], result$models[[2]])

}
