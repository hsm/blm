#' @export
blm <- function(model, prior = NULL, alpha = 1, beta = 1, ...) {
  frame <- model.frame(model, ...)
  m <- model.matrix(model, frame)

  if (is.null(prior)) {
    d <- dim(m)[[2]]
    prior <- list(mean = rep(0, d), covar = alpha * diag(1, nrow = d))
  }

  covar <- solve(prior$covar + beta * t(m) %*% m)
  mean <- beta * covar %*% t(m) %*% model.response(frame)
  posterior <- list(mean = mean, covar = covar)

  structure(list(formula = model,
                 frame = frame,
                 posterior = posterior,
                 call = match.call()),
            class = "blm")
}

distribution <- function(x) UseMethod("distribution")
distribution.default <- I
distribution.blm <- function(x) x$posterior

namevec <- function(x) {
  result <- as.vector(x)
  names(result) <- seq_along(x)
  result
}

coefficients <- function(object, ...) UseMethod("coefficients")

#' @export
coefficients.blm <- function(object, ...) {
  result <- as.vector(object$posterior$mean)
  names(result) <- rownames(object$posterior$mean)
  result
}

#' @export
predict.blm <- function(object, ...) {
  args <- list(...)
  d <- if (is.null(args$newdata)) object$frame else args$newdata
  predict_on_data(object, d)
}

predict_on_data <- function(object, data) {
  responseless.formula <- delete.response(terms(object$formula))
  frame <- model.frame(responseless.formula, data = data)
  m <- model.matrix(responseless.formula, frame)
  namevec(t(object$posterior$mean) %*% t(m))
}

#' @export
residuals.blm <- function(object, ...) {
  model.response(object$frame) - fitted(object, ...)
}

#' @export
fitted.blm <- function(object, ...) {
  predict_on_data(object, object$frame)
}

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

#' @export
deviance.blm <- function(object, ...) {
  sum(residuals(object, ...)^2)
}

#' @export
plot.blm <- function(object, ...) {
  r <- residuals(object, ...)
  yh <- fitted(object, ...)
  plot(yh, r, ylim = range(r), main = "Residuals vs Fitted", xlab = paste("Fitted", deparse(object$call), sep = "\n"), ylab = "Residuals")
  abline(h=0, col=4, lty=2)
}

#' @export
print.blm <- function(object, ...) {
  # Copied and modified from print.lm using getAnywhere('print.lm')
  cat("\nCall:\n", paste(deparse(object$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("Coefficients:\n")
  print.default(coefficients(object), print.gap = 2L, quote = FALSE)
  cat("\n")
  invisible(object)
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
