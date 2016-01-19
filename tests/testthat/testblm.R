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


