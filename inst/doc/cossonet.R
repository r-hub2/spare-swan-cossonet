## ----results = 'hide', message=FALSE, warning=FALSE---------------------------
devtools::install_github("jiieunshin/cossonet")
library(cossonet)
set.seed(20250101)

## -----------------------------------------------------------------------------
tr = data_generation(n = 200, p = 20, SNR = 9, response = "continuous")
str(tr)

te = data_generation(n = 1000, p = 20, SNR = 9, response = "continuous")
str(te)

## ----fig.width=8, fig.height=4------------------------------------------------
lambda0_seq = exp(seq(log(2^{-5}), log(2^{-1}), length.out = 20))
lambda_theta_seq = exp(seq(log(2^{-8}), log(2^{-5}), length.out = 20))

fit = cossonet(tr$x, tr$y, family = 'gaussian',
	       lambda0 = lambda0_seq,
	       lambda_theta = lambda_theta_seq
	       )

## -----------------------------------------------------------------------------
pred = cossonet.predict(fit, te$x)
str(pred)

mean((te$f - pred$f.new)^2)

