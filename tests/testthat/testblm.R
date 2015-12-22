context("Testing Bayesian linear regression")
data("test-data")

test_that("Prediction matches predicted values from test dataset", {
  expect_equal(length(x), 100)
  simple <- blm(y ~ x)
  prediction <- predict(simple)
  expect_equal(length(prediction), 100)
  expect_equal(prediction, predicted_y, tolerance = 1e-15)
})

test_that("Prediction works when data is added instead of being bound to formula", {
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



