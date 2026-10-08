# SER configurations (round 2: SER memo, table kernel, final SER step).
# Data: z x 6 for strong signals; R with a non-unit diagonal (a few
# distinct shat2); unstandardised XtX (no table); prior weights with zeros.
z6   <- d1$z * 6
Rd   <- d1$R; diag(Rd) <- 1 + (seq_len(nrow(Rd)) %% 3) * 1e-3
pw0  <- { w <- pw; w[seq(1, length(w), by = 10)] <- 0; w / sum(w) }
Xc   <- scale(d1$X, scale = FALSE)
ssu  <- list(XtX = crossprod(Xc), Xty = as.vector(crossprod(Xc, d1$y - mean(d1$y))),
             yty = sum((d1$y - mean(d1$y))^2), n = d1$n)
cases <- c(cases, list(
  ser_rss_z6       = function() susie_rss(z6, d1$R, n = d1$n),
  ser_rss_L20      = function() susie_rss(d1$z, d1$R, n = d1$n, L = 20),
  ser_rss_Rdiag    = function() susie_rss(d1$z, Rd, n = d1$n),
  ser_rss_pw0      = function() susie_rss(d1$z, d1$R, n = d1$n, prior_weights = pw0),
  ser_rss_nullw20  = function() susie_rss(d1$z, d1$R, n = d1$n, L = 20, null_weight = 0.5),
  ser_rss_zero     = function() susie_rss(rep(0, length(d1$z)), d1$R, n = d1$n),
  ser_rss_weak     = function() susie_rss(d1$z * 0.4, d1$R, n = d1$n, L = 5),
  ser_rss_zNA      = function() susie_rss(replace(d1$z, c(5, 50), NA), d1$R, n = d1$n),
  ser_rss_p1       = function() susie_rss(d1$z[1], d1$R[1, 1, drop = FALSE], n = d1$n, L = 2),
  ser_rss_V0       = function() susie_rss(d1$z, d1$R, n = d1$n, estimate_prior_variance = FALSE,
                                          scaled_prior_variance = 0),
  ser_rss_mism_mix = function() susie_rss(ds$z, ds$R, n = ds$n, R_mismatch = "eb_mix", R_finite = 500),
  ser_rss_mism150  = function() susie_rss(ds$z, cor(ds$X[1:150, ]), n = ds$n, R_mismatch = "eb", R_finite = 150),
  ser_rss_mism_wk  = function() susie_rss(ds$z * 0.3, ds$R, n = ds$n, L = 5, R_mismatch = "eb", R_finite = 300),
  ser_ss_unstd     = function() susie_ss(ssu$XtX, ssu$Xty, ssu$yty, ssu$n, standardize = FALSE),
  ser_ss_L20       = function() susie_ss(ss$XtX, ss$Xty, ss$yty, ss$n, L = 20),
  ser_ind_L20      = function() susie(d1$X, d1$y, L = 20),
  ser_ind_nullw    = function() susie(d1$X, d1$y, null_weight = 0.3, prior_weights = pw),
  ser_ind_strong   = function() susie(d1$X, d1$y + d1$X[, 200] * 0.6)
))
