library(testthat)
library(cossonet)

test_that("cossonet runs without errors", {
  set.seed(20250101)

  tr = data_generation(n = 50, p = 10, SNR = 9, response = "continuous")
  te = data_generation(n = 100, p = 10, SNR = 9, response = "continuous")

  lambda0_seq = exp(seq(log(2^-5), log(2^-1), length.out = 10))
  lambda_theta_seq = exp(seq(log(2^-8), log(2^-5), length.out = 10))

  fit = cossonet(tr$x, tr$y, family = 'gaussian',
                 lambda0 = lambda0_seq,
                 lambda_theta = lambda_theta_seq)

  pred = cossonet.predict(fit, te$x)

  expect_type(pred$f.new, "double")   # 예측값이 숫자형인지 확인
  expect_true(mean((te$f - pred$f.new)^2) >= 0)  # MSE가 음수가 아닌지 확인
})
