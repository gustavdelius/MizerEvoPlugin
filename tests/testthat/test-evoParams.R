test_that("evoParams works", {
  params <- evoParams()
  expect_is(params, "MizerParams")
})
