#' @export
blm <- function(model, prior = NULL, alpha = 1, beta = 1, ...) {
  frame <- model.frame(model, ...)

  if (is.null(prior)) {
    d <- dim(frame)[[2]]
    prior <- list(mean = rep(0, d), covar = alpha * diag(1, nrow = d))
  }

  # to update
  m <- model.matrix(model, frame)
  covar <- solve(prior$covar + beta * t(m) %*% m)
  mean <- beta * covar %*% t(m) %*% model.response(frame)
  posterior <- list(mean = mean, covar = covar)

  structure(list(formula = model,
                 frame = frame,
                 posterior = posterior),
            class = "blm")
}

distribution <- function(x) UseMethod("distribution")
distribution.default <- I
distribution.blm <- function(x) x$posterior

#update <- function(model, prior, ...) {
  #prior <- distribution(prior)
  # fitting code here
#}

namevec <- function(x) {
  result <- as.vector(x)
  names(result) <- seq_along(x)
  result
}

coefficients <- function(object, ...) UseMethod("coefficients")

#' @export
coefficients.blm <- function(object, ...) {
  namevec(object$posterior$mean)
}

#' @export
predict <- function(object, ...) UseMethod("predict")

#' @export
predict.blm <- function(object, ...) {
  responseless.formula <- delete.response(terms(object$formula))
  ellipsis_args <- list(...)
  frame <- model.frame(responseless.formula, data = ellipsis_args$newdata)
  m <- model.matrix(responseless.formula, frame)
  namevec(t(object$posterior$mean) %*% t(m))
}
