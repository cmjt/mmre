compile.mmre()
dll.mmre()
test_that("Simple two state continuous-time Markov model", {
  data(example)
  mod.basic <- fit.mmre(data = example$data,parameters = list(log_baseline = log(c(0.5,0.5))))
  p <- get.probs(mod.basic,1)
  expect_equal(round(p[1,1], 3), 0.790, tolerance = 0.01)
})
test_that("Model with individual level random effects", {
  data(example)
  mod.basic.re <- fit.mmre(data = example$data,parameters = example$parameters.basic.re)
  p <- get.probs(mod.basic.re,1)
  expect_equal(round(p[1,1], 3), 0.791, tolerance = 0.01)
})
test_that("Model with exponential decay effect of covariate and individual level random effects", {
  data(example)
  mod.decay.re <- fit.mmre(data = example$data,parameters = example$parameters.decay.re,
                           decay = TRUE, cov.names = "t.since")
  base <- get.coefs(mod.decay.re)$baseline_transition_matrix
  expect_equal(round(base[1,1], 3), -0.593, tolerance = 0.01)
})
