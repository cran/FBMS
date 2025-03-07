% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.R
\name{predict.bgnlm_model}
\alias{predict.bgnlm_model}
\title{Predict responses from a BGNLM model}
\usage{
\method{predict}{bgnlm_model}(
  object,
  x,
  link = function(x) {
     x
 },
  ...
)
}
\arguments{
\item{object}{A fitted \code{bgnlm_model} object obtained from the BGNLM fitting procedure.
It should contain the estimated coefficients in \code{model$coefs}.}

\item{x}{A \code{data.frame} containing the new data for which predictions are to be made.
The variables in \code{x} must match the features used in the model.}

\item{link}{A link function to apply to the linear predictor.
By default, it is the identity function \code{function(x)\{x\}},
but it can be any function such as \code{plogis} for logistic regression models.}

\item{...}{Additional arguments to pass to prediction function.}
}
\value{
A numeric vector of predicted values for the given data \code{x}.
These predictions are calculated as \eqn{\hat{y} = \text{link}(X \beta)},
where \eqn{X} is the design matrix and \eqn{\beta} are the model coefficients.
}
\description{
This function generates predictions from a fitted \code{bgnlm_model} object given a new dataset.
}
\examples{
\dontrun{
# Example with simulated data
set.seed(123)
x_train <- data.frame(PlanetaryMassJpt = rnorm(100), RadiusJpt = rnorm(100))
model <- list(
  coefs = c(Intercept = -0.5, PlanetaryMassJpt = 0.2, RadiusJpt = -0.1),
  class = "bgnlm_model"
)
class(model) <- "bgnlm_model"

# New data for prediction
x_new <- data.frame(PlanetaryMassJpt = c(0.1, -0.3), RadiusJpt = c(0.2, -0.1))

# Predict using the identity link (default)
preds <- predict.bgnlm_model(model, x_new)
}

}
