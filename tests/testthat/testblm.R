context("Testing Bayesian linear regression")
data("test-data")

test_that("Prediction matches predicted values from test dataset", {
  expect_equal(length(x), 100)
  simple <- blm(y ~ x)
  prediction <- predict(simple)
  expect_equal(length(prediction), 100)
  expect_equal(prediction, predicted_y, tolerance = 1e-15)
})

test_that("Prediction completes when formula does not have intercept", {
  simpler <- blm(y ~ x - 1)
  prediction <- predict(simpler)
  expect_equal(length(prediction), 100)
})

test_that("Prediction works when data is added instead of coming from global env", {
  df <- data.frame(a = y, b = x)
  simple <- blm(a ~ b, data = df)
  prediction <- predict(simple)
  expect_equal(length(prediction), 100)
  expect_equal(prediction, predicted_y, tolerance = 1e-15)
})

test_that("Updating in two steps ends roughly the same place", {
  df <- data.frame(a = y, b = x)
  simple <- blm(a ~ b, data = df)

  splitx <- split(x, ceiling(seq_along(x)/50))
  splity <- split(y, ceiling(seq_along(y)/50))
  df1 <- data.frame(a = splity$`1`, b = splitx$`1`)
  df2 <- data.frame(a = splity$`2`, b = splitx$`2`)
  first_half <- blm(a ~ b, data = df1)
  second_half <- update(first_half, data = df2)

  expect_equal(distribution(simple)$mean, distribution(second_half)$mean, tolerance=1e-1)
  expect_equal(distribution(simple)$covar, distribution(second_half)$covar, tolerance=1e-1)
  warning("One would expect the distributions to be more equal")
})

test_that("Prediction accepts new data", {
  simple <- blm(y ~ x)
  prediction <- predict(simple, newdata = data.frame(x = c(0.59, 0.876)))
  expect_equal(length(prediction), 2)
  expect_equal(as.vector(prediction), c(1.516,1.883), tolerance=1e-3)
})

test_that("Fitted does not accept new data", {
  simple <- blm(y ~ x)
  fit <- fitted(simple, newdata = data.frame(x = c(0.59, 0.876)))
  expect_equal(length(fit), 100)
})

test_that("Deviance equals the sum of square distances", {
  simple <- blm(y ~ x)
  expect_equal(deviance(simple), 8.62412, tolerance=1e-5)
})

test_that("Confidence intervals are symmetric around the means", {
  simple <- blm(y ~ x)
  conf <- confint(simple, level = 0.90)
  expect_equal(rowMeans(conf), coefficients(simple))
})

make_prior <- function(alpha = 1) {
  list(mu = c(0,0), sigma = alpha * diag(1, nrow = 2))
}

sample_from_prior <- function(n, alpha = 1) {
  prior <- make_prior(alpha)
  df <- as.data.frame(mvrnorm(n = n, mu = prior$mu, Sigma = prior$sigma))
  colnames(df) <- c("w0", "w1")
  df
}

make_posterior <- function(x, y, alpha = 1, beta = 1) {
  Phi = cbind(1, x)
  sxy = solve(alpha * diag(1, nrow = 2) + beta * t(Phi) %*% Phi)
  mxy = beta * sxy %*% t(Phi) %*% y
  list(mu = as.vector(mxy), sigma = sxy)
}

sample_from_posterior <- function(n, x, y, alpha = 1, beta = 1) {
  posterior <- make_posterior(x, y, alpha, beta)
  df <- as.data.frame(mvrnorm(n = n, mu = posterior$mu, Sigma = posterior$sigma))
  colnames(df) <- c("w0", "w1")
  df
}

sse <- function(true_weights, derived_weights) {
  sum((true_weights - derived_weights)^2)
}

seed <- as.integer(1000 * rnorm(1))
test_that(paste("more data improves the model with seed", seed), {
  set.seed(seed)
  w0 <- 0.3 ; w1 <- 1.1 ; beta <- 1

  x <- rnorm(10)
  y <- rnorm(10, w1 * x + w0, 1/beta)
  obj <- blm(y ~ x)

  current_error <- sse(c(w0, w1), distribution(obj)$mean)
  current_var_magnitude <- sum(diag((distribution(obj)$covar)^2))

  x <- rnorm(1000)
  y <- rnorm(1000, w1 * x + w0, 1/beta)
  obj <- blm(y ~ x)

  prev_error <- current_error
  current_error <- sse(c(w0, w1), distribution(obj)$mean)
  expect_true(current_error < prev_error)

  prev_var_magnitude <- current_var_magnitude
  current_var_magnitude <- sum(diag((distribution(obj)$covar)^2))
  expect_true(current_var_magnitude < prev_var_magnitude)

  x <- rnorm(10000)
  y <- rnorm(10000, w1 * x + w0, 1/beta)
  obj <- blm(y ~ x)

  prev_error <- current_error
  current_error <- sse(c(w0, w1), distribution(obj)$mean)
  expect_true(current_error < prev_error)

  prev_var_magnitude <- current_var_magnitude
  current_var_magnitude <- sum(diag((distribution(obj)$covar)^2))
  expect_true(current_var_magnitude < prev_var_magnitude)
})
