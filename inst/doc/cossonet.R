## ----eval=FALSE---------------------------------------------------------------
# install.packages("cossonet")
# library(cossonet)
# set.seed(20250101)

## ----eval=FALSE---------------------------------------------------------------
# tr = data_generation(n = 200, p = 20, SNR = 9, response = "continuous")
# te = data_generation(n = 1000, p = 20, SNR = 9, response = "continuous")

## ----eval=FALSE---------------------------------------------------------------
# lambda0_seq = exp(seq(log(2^{-5}), log(2^{-1}), length.out = 20))
# lambda_theta_seq = exp(seq(log(2^{-8}), log(2^{-5}), length.out = 20))
# 
# fit = cossonet(tr$x, tr$y, family = 'gaussian',
# 	       lambda0 = lambda0_seq,
# 	       lambda_theta = lambda_theta_seq
# 	       )

## ----eval=FALSE---------------------------------------------------------------
# pred = cossonet.predict(fit, te$x)
# mean((te$f - pred$f.new)^2)

