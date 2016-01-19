#' Fitting Bayesian linear models
#'
#' @param formula An object of class \link{formula}
#' @param prior a prior
#' @param alpha the covariance precision
#' @param beta the precision
#' @param ... other arguments passed to the model.frame method. Adding data as a parameter can be used to add data not in the enclosing environment
#' @return an object of \link{class} blm
#' @export
blm <- function(formula, prior = NULL, alpha = 1, beta = 1, ...) {
  frame <- model.frame(formula, ...)
  m <- model.matrix(formula, frame)

  if (is.null(prior)) {
    d <- dim(m)[[2]]
    prior <- list(mean = rep(0, d), covar = alpha * diag(1, nrow = d))
  }

  covar <- solve(prior$covar + beta * t(m) %*% m)
  mean <- beta * covar %*% t(m) %*% model.response(frame)
  posterior <- list(mean = mean, covar = covar)

  structure(list(formula = formula,
                 frame = frame,
                 posterior = posterior,
                 call = sys.call()),
            class = "blm")
}

distribution <- function(object) UseMethod("distribution")
distribution.default <- I

#' Extract model distribution
#'
#' The distribution parameters is the means and the covariance of the current posterior
#'
#' @param object an object of \link{class} blm
#' @return The distribution parameters
#' @export
distribution.blm <- function(object) object$posterior

update <- function(object, ...) UseMethod("update")
update.default <- function(object, ...) stop("update default not implemented")

#' Updating the model
#'
#' The distribution of the current model is used as a prior, and any new data is used to build the
#' new model.
#'
#' @param object the current model. An object of \link{class} blm
#' @param ... other arguments passed to blm
#' @return An updated object of \link{class} blm
#' @export
update.blm <- function(object, ...) {
  blm(object$formula, prior = distribution(object), ...)
}

coefficients <- function(object, ...) UseMethod("coefficients")

#' Extract model coefficients
#'
#' The coefficients of the blm model is the means of the posterior
#'
#' @param object an object of \link{class} blm
#' @param ... other arguments
#' @return A named numeric vector of coefficients
#' @export
coefficients.blm <- function(object, ...) {
  result <- as.vector(object$posterior$mean)
  names(result) <- rownames(object$posterior$mean)
  result
}

#' Predict method for Bayesian Linear models
#'
#' Provides prediction based on the fitted Bayesian linear model
#'
#' @param object an object of \link{class} blm
#' @param newdata An optional data frame in which to look for variables with which to predict. If omitted, the fitted values are used
#' @param ... other arguments
#' @return A vector of predictions
#' @export
predict.blm <- function(object, newdata = NULL, ...) {
  d <- if (is.null(newdata)) object$frame else newdata
  predict_on_data(object, d)
}

predict_on_data <- function(object, data) {
  responseless.formula <- delete.response(terms(object$formula))
  frame <- model.frame(responseless.formula, data = data)
  m <- model.matrix(responseless.formula, frame)
  namevec(t(object$posterior$mean) %*% t(m))
}

namevec <- function(x) {
  result <- as.vector(x)
  names(result) <- seq_along(x)
  result
}

#' Extract blm residuals
#'
#' Calculates the difference between the predicted values and the observed values for the response variable
#'
#' @param object an object of \link{class} blm
#' @param ... other arguments
#' @return A vector of residuals
#' @export
residuals.blm <- function(object, ...) {
  model.response(object$frame) - fitted(object, ...)
}

#' Extract blm Fitted Values
#'
#' The fitted values are the predicted values of the model using the data used to fit the model
#'
#' @param object an object of \link{class} blm
#' @param ... other arguments
#' @return A vector of fits
#' @export
fitted.blm <- function(object, ...) {
  predict_on_data(object, object$frame)
}

#' Confidence Intervals for Model Parameters
#'
#' Computes confidence intervals for the parameters in the fitted model
#'
#' @param object an object of \link{class} blm
#' @param level the confidence level required
#' @param ... other arguments
#' @return A matrix with columns giving lower and upper confidence limits for each parameter.
#' @export
confint.blm <- function(object, level = 0.95, ...) {
  alpha <- 1 - level
  mean <- object$posterior$mean
  sd <- sqrt(diag(object$posterior$covar))
  critval <- qnorm(1 - alpha / 2, mean, sd)
  lowerconf <- mean - critval * sd
  upperconf <- mean + critval * sd
  result <- cbind(lowerconf, upperconf)
  colnames(result) <- c(paste(100 * alpha / 2, "%"), paste(100 * (1 - alpha / 2), "%"))
  result
}

#' Model Deviance
#'
#' Computes the sum of the squared distance between the observed and predicted values
#'
#' @param object an object of \link{class} blm
#' @param ... other arguments
#' @return the deviance
#' @export
deviance.blm <- function(object, ...) {
  sum(residuals(object, ...)^2)
}

#' Plot Diagnostics for a blm Object
#'
#' Plots residuals vs fitted parameters
#'
#' @param x an object of \link{class} blm
#' @param ... other arguments
#' @export
plot.blm <- function(x, ...) {
  r <- residuals(x, ...)
  yh <- fitted(x, ...)
  plot(yh, r, ylim = range(r), main = "Residuals vs Fitted", xlab = paste("Fitted", deparse(x$call), sep = "\n"), ylab = "Residuals")
  abline(h=0, col=4, lty=2)
}

#' @export
print.blm <- function(x, ...) {
  # Copied and modified from print.lm using getAnywhere('print.lm')
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("Coefficients:\n")
  print.default(coefficients(x), print.gap = 2L, quote = FALSE)
  cat("\n")
  invisible(x)
}

#' @export
summary.blm <- function(object, ...) {
  print.blm(object, ...)
  cat("Covariance matrix:\n")
  print.default(object$posterior$covar, print.gap = 2L, quote = FALSE)
  cat("\n")
  cat("Deviance: ", deviance(object, ...), "\n\n")
  cat("Confidence at 95% level:\n")
  print.default(confint(object), print.gap = 2L, quote = FALSE)
  cat("\n")
}
