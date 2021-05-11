library(fmsb)

context("Basic statistics")

test_that("SIQR is semi-interquartile range", {
	expect_equal(SIQR(c(1, 5, 10, 100, 1000)), 47.5)
	expect_equal(SIQR(1:100), 25)
})

test_that("True median", {
  expect_equal(truemedian(c(5,6,6,7,7,8,8,8)), 7)
  expect_equal(truemedian(c(3,3,4,4,4)), 11/3)
})

test_that("Geary's test for normality", {
  expect_equal(geary.test(20:50)$p.value, 0.96207213)
})

context("Epidemiologic estimates")

test_that("ratedifference", {
  expect_equal(ratedifference(136, 1709, 22050, 127650, CRC=TRUE)$estimate, -0.00722037)
})