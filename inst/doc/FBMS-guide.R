## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(FBMS)
library(GenSA)
library(fastglm)
data("breastcancer")
bc <- breastcancer[,c(ncol(breastcancer),2:(ncol(breastcancer)-1))]

## -----------------------------------------------------------------------------
to3 <- function(x) x^3

transforms <- c("sigmoid","sin_deg","exp_dbl","p0","troot","to3")

## -----------------------------------------------------------------------------
probs <- gen.probs.gmjmcmc(transforms)
params <- gen.params.gmjmcmc(bc)

## -----------------------------------------------------------------------------
loglik.example <- function (y, x, model, complex, params) {
  r <- 20/223
  suppressWarnings({mod <- fastglm(as.matrix(x[,model]), y, family=binomial())})
  ret <- (-(mod$deviance -2*log(r)*sum(complex$width)))/2
  return(list(crit=ret, coefs=mod$coefficients))
}

## -----------------------------------------------------------------------------
params$feat$alpha <- 3

## -----------------------------------------------------------------------------
set.seed(1234)
result <- gmjmcmc(bc, loglik.example, logistic.loglik.alpha, transforms, P=3, N=30, N.final=60, probs, params)

## ---- fig.width=6, fig.height=6-----------------------------------------------
plot(result)

