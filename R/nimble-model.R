boundedSliceSampler <- nimbleFunction(
  contains = nimble::sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    # Create internal slice sampler safely
    sliceSampler <- nimble::sampler_slice(
      model = model,
      mvSaved = mvSaved,
      target = target,
      control = control
    )

    # Optionally accept user-supplied K from control
    if (!is.null(control$K)) {
      K <- control$K
    } else {
      # Default: try to detect length of p if it exists
      if ("p" %in% model$getNodeNames()) {
        K <- length(model[['p']])
      } else {
        stop(
          "Please provide K (number of categories) in 'control' for boundedSliceSampler."
        )
      }
    }
  },
  run = function() {
    sliceSampler$run()

    value <- model[[target]]

    # Apply integer bounds
    if (value < 1) {
      value <- 1
    }
    if (value > K) {
      value <- K
    }

    model[[target]] <<- value
  },
  methods = list(
    reset = function() {
      sliceSampler$reset()
    }
  )
)


nimblemod <- nimbleCode({
  col.ind ~ dcat(w[1:K])
  for (k in 1:K) {
    vec[k] <- equals(col.ind, k)
  }

  if (family != "multinomial") {
    if (!(family %in% c("prophaz", "clogit"))) {
      a0 ~ dnorm(0, 0.001)
    } else {
      a0 <- 0
    }

    if (Xlen > 1) {
      for (k in 1:Xlen) {
        a[k] ~ dnorm(0, 0.001)
      }
    } else if (Xlen == 1) {
      a ~ dnorm(0, .001)
    }

    if (doseRRmod == "EXP" | family == "gaussian") {
      # Normal priors with large variance
      if (deg > 1) {
        for (k in 1:deg) {
          b[k] ~ dnorm(0, 0.001)
        }
      } else {
        b ~ dnorm(0, .001)
      }

      if (Mlen > 1) {
        if (deg == 1) {
          for (k in 1:Mlen) {
            bm[k] ~ dnorm(0, 0.001)
          }
        } else if (deg > 1) {
          for (k in 1:Mlen) {
            bm1[k] ~ dnorm(0, 0.001)
            bm2[k] ~ dnorm(0, 0.001)
          }
        }
      } else if (Mlen == 1) {
        if (deg == 1) {
          bm ~ dnorm(0, 0.001)
        } else if (deg > 1) {
          bm1 ~ dnorm(0, 0.001)
          bm2 ~ dnorm(0, 0.001)
        }
      }
    } else if (doseRRmod == "ERR") {
      if (ERRprior == "truncated_horseshoe") {
        tau ~ T(dt(0, 1, df = 1), 0, )
        if (deg > 1) {
          for (kk in 1:2) {
            # horseshoe prior
            lambda[kk] ~ T(dt(0, 1, df = 1), 0, )
            b[kk] ~ T(dnorm(0, sd = lambda[kk] * tau), 0, )
          }
        } else {
          lambda ~ T(dt(0, 1, df = 1), 0, )
          b ~ T(dnorm(0, sd = lambda * tau), 0, )
        }

        if (Mlen > 1) {
          if (deg == 1) {
            for (zz in 1:Mlen) {
              lambda_m[zz] ~ T(dt(0, 1, df = 1), 0, )
              bm[zz] ~ T(dnorm(0, sd = lambda_m[zz] * tau), 0, )
            }
          } else if (deg > 1) {
            for (zz in 1:Mlen) {
              # horseshoe prior
              lambda_m1[zz] ~ T(dt(0, 1, df = 1), 0, )
              bm1[zz] ~ T(dnorm(0, sd = lambda_m1[zz] * tau), 0, )
              lambda_m2[zz] ~ T(dt(0, 1, df = 1), 0, )
              bm2[zz] ~ T(dnorm(0, sd = lambda_m2[zz] * tau), 0, )
            }
          }
        } else if (Mlen == 1) {
          if (deg == 1) {
            lambda_m ~ T(dt(0, 1, df = 1), 0, )
            bm ~ T(dnorm(0, sd = lambda_m * tau), 0, )
          } else if (deg > 1) {
            lambda_m1 ~ T(dt(0, 1, df = 1), 0, )
            bm1 ~ T(dnorm(0, sd = lambda_m1 * tau), 0, )
            lambda_m2 ~ T(dt(0, 1, df = 1), 0, )
            bm2 ~ T(dnorm(0, sd = lambda_m2 * tau), 0, )
          }
        }
      } else if (ERRprior == "truncated_normal") {
        if (deg > 1) {
          for (k in 1:deg) {
            b[k] ~ T(dnorm(0, 0.001), 0, )
          }
        } else {
          b ~ T(dnorm(0, .001), 0, )
        }

        if (Mlen > 1) {
          if (deg == 1) {
            for (k in 1:Mlen) {
              bm[k] ~ T(dnorm(0, 0.001), 0, )
            }
          } else if (deg > 1) {
            for (k in 1:Mlen) {
              bm1[k] ~ T(dnorm(0, 0.001), 0, )
              bm2[k] ~ T(dnorm(0, 0.001), 0, )
            }
          }
        } else if (Mlen == 1) {
          if (deg == 1) {
            bm ~ T(dnorm(0, 0.001), 0, )
          } else if (deg > 1) {
            bm1 ~ T(dnorm(0, 0.001), 0, )
            bm2 ~ T(dnorm(0, 0.001), 0, )
          }
        }
      } else if (ERRprior == "truncated_doubleexponential") {
        if (deg > 1) {
          for (kk in 1:2) {
            lambda[kk] ~ T(dt(0, 1, df = 1), 0, )
            b[kk] ~ T(ddexp(0.0, lambda[kk]), 0, )
          }
        } else {
          lambda ~ T(dt(0, 1, df = 1), 0, )
          b ~ T(ddexp(0.0, lambda), 0, )
        }

        if (Mlen > 1) {
          if (deg == 1) {
            for (zz in 1:Mlen) {
              lambda_m[zz] ~ T(dt(0, 1, df = 1), 0, )
              bm[zz] ~ T(ddexp(0.0, lambda_m[zz]), 0, )
            }
          } else if (deg > 1) {
            for (zz in 1:Mlen) {
              lambda_m1[zz] ~ T(dt(0, 1, df = 1), 0, )
              bm1[zz] ~ T(ddexp(0.0, lambda_m1[zz]), 0, )
              lambda_m2[zz] ~ T(dt(0, 1, df = 1), 0, )
              bm2[zz] ~ T(ddexp(0.0, lambda_m2[zz]), 0, )
            }
          }
        } else if (Mlen == 1) {
          if (deg == 1) {
            lambda_m ~ T(dt(0, 1, df = 1), 0, )
            bm ~ T(ddexp(0.0, lambda_m), 0, )
          } else if (deg > 1) {
            lambda_m1 ~ T(dt(0, 1, df = 1), 0, )
            bm1 ~ T(ddexp(0.0, lambda_m1), 0, )
            lambda_m2 ~ T(dt(0, 1, df = 1), 0, )
            bm2 ~ T(ddexp(0.0, lambda_m2), 0, )
          }
        }
      } else if (ERRprior == "horseshoe") {
        tau ~ T(dt(0, 1, df = 1), 0, )
        if (deg > 1) {
          for (kk in 1:2) {
            # horseshoe prior
            lambda[kk] ~ T(dt(0, 1, df = 1), 0, )
            b[kk] ~ dnorm(0, sd = lambda[kk] * tau)
          }
        } else {
          lambda ~ T(dt(0, 1, df = 1), 0, )
          b ~ dnorm(0, sd = lambda * tau)
        }

        if (Mlen > 1) {
          if (deg == 1) {
            for (zz in 1:Mlen) {
              lambda_m[zz] ~ T(dt(0, 1, df = 1), 0, )
              bm[zz] ~ dnorm(0, sd = lambda_m[zz] * tau)
            }
          } else if (deg > 1) {
            for (zz in 1:Mlen) {
              # horseshoe prior
              lambda_m1[zz] ~ T(dt(0, 1, df = 1), 0, )
              bm1[zz] ~ dnorm(0, sd = lambda_m1[zz] * tau)
              lambda_m2[zz] ~ T(dt(0, 1, df = 1), 0, )
              bm2[zz] ~ dnorm(0, sd = lambda_m2[zz] * tau)
            }
          }
        } else if (Mlen == 1) {
          if (deg == 1) {
            lambda_m ~ T(dt(0, 1, df = 1), 0, )
            bm ~ dnorm(0, sd = lambda_m * tau)
          } else if (deg > 1) {
            lambda_m1 ~ T(dt(0, 1, df = 1), 0, )
            bm1 ~ dnorm(0, sd = lambda_m1 * tau)
            lambda_m2 ~ T(dt(0, 1, df = 1), 0, )
            bm2 ~ dnorm(0, sd = lambda_m2 * tau)
          }
        }
      } else if (ERRprior == "normal") {
        if (deg > 1) {
          for (k in 1:deg) {
            b[k] ~ dnorm(0, 0.001)
          }
        } else {
          b ~ dnorm(0, .001)
        }

        if (Mlen > 1) {
          if (deg == 1) {
            for (k in 1:Mlen) {
              bm[k] ~ dnorm(0, 0.001)
            }
          } else if (deg > 1) {
            for (k in 1:Mlen) {
              bm1[k] ~ dnorm(0, 0.001)
              bm2[k] ~ dnorm(0, 0.001)
            }
          }
        } else if (Mlen == 1) {
          if (deg == 1) {
            bm ~ dnorm(0, 0.001)
          } else if (deg > 1) {
            bm1 ~ dnorm(0, 0.001)
            bm2 ~ dnorm(0, 0.001)
          }
        }
      } else if (ERRprior == "doubleexponential") {
        if (deg > 1) {
          for (kk in 1:2) {
            lambda[kk] ~ T(dt(0, 1, df = 1), 0, )
            b[kk] ~ ddexp(0.0, lambda[kk])
          }
        } else {
          lambda ~ T(dt(0, 1, df = 1), 0, )
          b ~ ddexp(0.0, lambda)
        }

        if (Mlen > 1) {
          if (deg == 1) {
            for (zz in 1:Mlen) {
              lambda_m[zz] ~ T(dt(0, 1, df = 1), 0, )
              bm[zz] ~ ddexp(0.0, lambda_m[zz])
            }
          } else if (deg > 1) {
            for (zz in 1:Mlen) {
              lambda_m1[zz] ~ T(dt(0, 1, df = 1), 0, )
              bm1[zz] ~ ddexp(0.0, lambda_m1[zz])
              lambda_m2[zz] ~ T(dt(0, 1, df = 1), 0, )
              bm2[zz] ~ ddexp(0.0, lambda_m2[zz])
            }
          }
        } else if (Mlen == 1) {
          if (deg == 1) {
            lambda_m ~ T(dt(0, 1, df = 1), 0, )
            bm ~ ddexp(0.0, lambda_m)
          } else if (deg > 1) {
            lambda_m1 ~ T(dt(0, 1, df = 1), 0, )
            bm1 ~ ddexp(0.0, lambda_m1)
            lambda_m2 ~ T(dt(0, 1, df = 1), 0, )
            bm2 ~ ddexp(0.0, lambda_m2)
          }
        }
      }
    } else if (doseRRmod == "LINEXP") {
      if (ERRprior == "truncated_horseshoe") {
        tau ~ T(dt(0, 1, df = 1), 0, )

        # horseshoe prior
        lambda[1] ~ T(dt(0, 1, df = 1), 0, )
        b[1] ~ T(dnorm(0, sd = lambda[1] * tau), 0, )
        lambda[2] ~ T(dt(0, 1, df = 1), 0, )
        b[2] ~ dnorm(0, sd = lambda[2] * tau)

        if (Mlen > 1) {
          for (zz in 1:Mlen) {
            # horseshoe prior
            lambda_m1[zz] ~ T(dt(0, 1, df = 1), 0, )
            bm1[zz] ~ T(dnorm(0, sd = lambda_m1[zz] * tau), 0, )
            lambda_m2[zz] ~ T(dt(0, 1, df = 1), 0, )
            bm2[zz] ~ dnorm(0, sd = lambda_m2[zz] * tau)
          }
        } else if (Mlen == 1) {
          lambda_m1 ~ T(dt(0, 1, df = 1), 0, )
          bm1 ~ T(dnorm(0, sd = lambda_m1 * tau), 0, )
          lambda_m2 ~ T(dt(0, 1, df = 1), 0, )
          bm2 ~ dnorm(0, sd = lambda_m2 * tau)
        }
      } else if (ERRprior == "truncated_normal") {
        b[1] ~ T(dnorm(0, 0.001), 0, )
        b[2] ~ dnorm(0, 0.001)

        if (Mlen > 1) {
          for (k in 1:Mlen) {
            bm1[k] ~ T(dnorm(0, 0.001), 0, )
            bm2[k] ~ dnorm(0, 0.001)
          }
        } else if (Mlen == 1) {
          bm1 ~ T(dnorm(0, 0.001), 0, )
          bm2 ~ dnorm(0, 0.001)
        }
      } else if (ERRprior == "truncated_doubleexponential") {
        lambda[1] ~ T(dt(0, 1, df = 1), 0, )
        b[1] ~ T(ddexp(0.0, lambda[1]), 0, )
        lambda[2] ~ T(dt(0, 1, df = 1), 0, )
        b[2] ~ ddexp(0.0, lambda[2])

        if (Mlen > 1) {
          for (zz in 1:Mlen) {
            lambda_m1[zz] ~ T(dt(0, 1, df = 1), 0, )
            bm1[zz] ~ T(ddexp(0.0, lambda_m1[zz]), 0, )
            lambda_m2[zz] ~ T(dt(0, 1, df = 1), 0, )
            bm2[zz] ~ ddexp(0.0, lambda_m2[zz])
          }
        } else if (Mlen == 1) {
          lambda_m1 ~ T(dt(0, 1, df = 1), 0, )
          bm1 ~ T(ddexp(0.0, lambda_m1), 0, )
          lambda_m2 ~ T(dt(0, 1, df = 1), 0, )
          bm2 ~ ddexp(0.0, lambda_m2)
        }
      } else if (ERRprior == "horseshoe") {
        tau ~ T(dt(0, 1, df = 1), 0, )

        for (kk in 1:2) {
          # horseshoe prior
          lambda[kk] ~ T(dt(0, 1, df = 1), 0, )
          b[kk] ~ dnorm(0, sd = lambda[kk] * tau)
        }

        if (Mlen > 1) {
          for (zz in 1:Mlen) {
            # horseshoe prior
            lambda_m1[zz] ~ T(dt(0, 1, df = 1), 0, )
            bm1[zz] ~ dnorm(0, sd = lambda_m1[zz] * tau)
            lambda_m2[zz] ~ T(dt(0, 1, df = 1), 0, )
            bm2[zz] ~ dnorm(0, sd = lambda_m2[zz] * tau)
          }
        } else if (Mlen == 1) {
          lambda_m1 ~ T(dt(0, 1, df = 1), 0, )
          bm1 ~ dnorm(0, sd = lambda_m1 * tau)
          lambda_m2 ~ T(dt(0, 1, df = 1), 0, )
          bm2 ~ dnorm(0, sd = lambda_m2 * tau)
        }
      } else if (ERRprior == "normal") {
        for (k in 1:2) {
          b[k] ~ dnorm(0, 0.001)
        }

        if (Mlen > 1) {
          for (k in 1:Mlen) {
            bm1[k] ~ dnorm(0, 0.001)
            bm2[k] ~ dnorm(0, 0.001)
          }
        } else if (Mlen == 1) {
          bm1 ~ dnorm(0, 0.001)
          bm2 ~ dnorm(0, 0.001)
        }
      } else if (ERRprior == "doubleexponential") {
        for (kk in 1:2) {
          lambda[kk] ~ T(dt(0, 1, df = 1), 0, )
          b[kk] ~ ddexp(0.0, lambda[kk])
        }

        if (Mlen > 1) {
          for (zz in 1:Mlen) {
            lambda_m1[zz] ~ T(dt(0, 1, df = 1), 0, )
            bm1[zz] ~ ddexp(0.0, lambda_m1[zz])
            lambda_m2[zz] ~ T(dt(0, 1, df = 1), 0, )
            bm2[zz] ~ ddexp(0.0, lambda_m2[zz])
          }
        } else if (Mlen == 1) {
          lambda_m1 ~ T(dt(0, 1, df = 1), 0, )
          bm1 ~ ddexp(0.0, lambda_m1)
          lambda_m2 ~ T(dt(0, 1, df = 1), 0, )
          bm2 ~ ddexp(0.0, lambda_m2)
        }
      }
    }
  } else {
    # Multinomial

    for (z in 1:(Z - 1)) {
      a0[z] ~ dnorm(0, 0.001)

      if (Xlen > 1) {
        for (k in 1:Xlen) {
          a[k, z] ~ dnorm(0, 0.001)
        }
      } else if (Xlen == 1) {
        a[z] ~ dnorm(0, .001)
      }
    }

    if (doseRRmod == "EXP") {
      # Normal priors with large variance
      if (deg > 1) {
        for (z in 1:(Z - 1)) {
          for (k in 1:deg) {
            b[k, z] ~ dnorm(0, 0.001)
          }
        }
      } else {
        for (z in 1:(Z - 1)) {
          b[z] ~ dnorm(0, .001)
        }
      }

      if (Mlen > 1) {
        if (deg == 1) {
          for (z in 1:(Z - 1)) {
            for (k in 1:Mlen) {
              bm[k, z] ~ dnorm(0, 0.001)
            }
          }
        } else if (deg > 1) {
          for (z in 1:(Z - 1)) {
            for (k in 1:Mlen) {
              bm1[k, z] ~ dnorm(0, 0.001)
              bm2[k, z] ~ dnorm(0, 0.001)
            }
          }
        }
      } else if (Mlen == 1) {
        if (deg == 1) {
          for (z in 1:(Z - 1)) {
            bm[z] ~ dnorm(0, 0.001)
          }
        } else if (deg > 1) {
          for (z in 1:(Z - 1)) {
            bm1[z] ~ dnorm(0, 0.001)
            bm2[z] ~ dnorm(0, 0.001)
          }
        }
      }
    } else if (doseRRmod == "ERR") {
      if (ERRprior == "truncated_horseshoe") {
        tau ~ T(dt(0, 1, df = 1), 0, )
        if (deg > 1) {
          for (z in 1:(Z - 1)) {
            for (kk in 1:2) {
              # horseshoe prior
              lambda[kk, z] ~ T(dt(0, 1, df = 1), 0, )
              b[kk, z] ~ T(dnorm(0, sd = lambda[kk, z] * tau), 0, )
            }
          }
        } else {
          for (z in 1:(Z - 1)) {
            lambda[z] ~ T(dt(0, 1, df = 1), 0, )
            b[z] ~ T(dnorm(0, sd = lambda[z] * tau), 0, )
          }
        }

        if (Mlen > 1) {
          if (deg == 1) {
            for (z in 1:(Z - 1)) {
              for (zz in 1:Mlen) {
                lambda_m[zz, z] ~ T(dt(0, 1, df = 1), 0, )
                bm[zz, z] ~ T(dnorm(0, sd = lambda_m[zz, z] * tau), 0, )
              }
            }
          } else if (deg > 1) {
            for (z in 1:(Z - 1)) {
              for (zz in 1:Mlen) {
                # horseshoe prior
                lambda_m1[zz, z] ~ T(dt(0, 1, df = 1), 0, )
                bm1[zz, z] ~ T(dnorm(0, sd = lambda_m1[zz, z] * tau), 0, )
                lambda_m2[zz, z] ~ T(dt(0, 1, df = 1), 0, )
                bm2[zz, z] ~ T(dnorm(0, sd = lambda_m2[zz, z] * tau), 0, )
              }
            }
          }
        } else if (Mlen == 1) {
          if (deg == 1) {
            for (z in 1:(Z - 1)) {
              lambda_m[z] ~ T(dt(0, 1, df = 1), 0, )
              bm[z] ~ T(dnorm(0, sd = lambda_m[z] * tau), 0, )
            }
          } else if (deg > 1) {
            for (z in 1:(Z - 1)) {
              lambda_m1[z] ~ T(dt(0, 1, df = 1), 0, )
              bm1[z] ~ T(dnorm(0, sd = lambda_m1[z] * tau), 0, )
              lambda_m2[z] ~ T(dt(0, 1, df = 1), 0, )
              bm2[z] ~ T(dnorm(0, sd = lambda_m2[z] * tau), 0, )
            }
          }
        }
      } else if (ERRprior == "truncated_normal") {
        if (deg > 1) {
          for (z in 1:(Z - 1)) {
            for (k in 1:deg) {
              b[k, z] ~ T(dnorm(0, 0.001), 0, )
            }
          }
        } else {
          for (z in 1:(Z - 1)) {
            b[z] ~ T(dnorm(0, .001), 0, )
          }
        }

        if (Mlen > 1) {
          if (deg == 1) {
            for (z in 1:(Z - 1)) {
              for (k in 1:Mlen) {
                bm[k, z] ~ T(dnorm(0, 0.001), 0, )
              }
            }
          } else if (deg > 1) {
            for (z in 1:(Z - 1)) {
              for (k in 1:Mlen) {
                bm1[k, z] ~ T(dnorm(0, 0.001), 0, )
                bm2[k, z] ~ T(dnorm(0, 0.001), 0, )
              }
            }
          }
        } else if (Mlen == 1) {
          if (deg == 1) {
            for (z in 1:(Z - 1)) {
              bm[z] ~ T(dnorm(0, 0.001), 0, )
            }
          } else if (deg > 1) {
            for (z in 1:(Z - 1)) {
              bm1[z] ~ T(dnorm(0, 0.001), 0, )
              bm2[z] ~ T(dnorm(0, 0.001), 0, )
            }
          }
        }
      } else if (ERRprior == "truncated_doubleexponential") {
        if (deg > 1) {
          for (z in 1:(Z - 1)) {
            for (kk in 1:2) {
              lambda[kk, z] ~ T(dt(0, 1, df = 1), 0, )
              b[kk, z] ~ T(ddexp(0.0, lambda[kk, z]), 0, )
            }
          }
        } else {
          for (z in 1:(Z - 1)) {
            lambda[z] ~ T(dt(0, 1, df = 1), 0, )
            b[z] ~ T(ddexp(0.0, lambda), 0, )
          }
        }

        if (Mlen > 1) {
          if (deg == 1) {
            for (z in 1:(Z - 1)) {
              for (zz in 1:Mlen) {
                lambda_m[zz, z] ~ T(dt(0, 1, df = 1), 0, )
                bm[zz, z] ~ T(ddexp(0.0, lambda_m[zz, z]), 0, )
              }
            }
          } else if (deg > 1) {
            for (z in 1:(Z - 1)) {
              for (zz in 1:Mlen) {
                lambda_m1[zz, z] ~ T(dt(0, 1, df = 1), 0, )
                bm1[zz, z] ~ T(ddexp(0.0, lambda_m1[zz, z]), 0, )
                lambda_m2[zz, z] ~ T(dt(0, 1, df = 1), 0, )
                bm2[zz, z] ~ T(ddexp(0.0, lambda_m2[zz, z]), 0, )
              }
            }
          }
        } else if (Mlen == 1) {
          if (deg == 1) {
            for (z in 1:(Z - 1)) {
              lambda_m[z] ~ T(dt(0, 1, df = 1), 0, )
              bm[z] ~ T(ddexp(0.0, lambda_m[z]), 0, )
            }
          } else if (deg > 1) {
            for (z in 1:(Z - 1)) {
              lambda_m1[z] ~ T(dt(0, 1, df = 1), 0, )
              bm1[z] ~ T(ddexp(0.0, lambda_m1[z]), 0, )
              lambda_m2[z] ~ T(dt(0, 1, df = 1), 0, )
              bm2[z] ~ T(ddexp(0.0, lambda_m2[z]), 0, )
            }
          }
        }
      } else if (ERRprior == "horseshoe") {
        tau ~ T(dt(0, 1, df = 1), 0, )
        if (deg > 1) {
          for (z in 1:(Z - 1)) {
            for (kk in 1:2) {
              # horseshoe prior
              lambda[kk, z] ~ T(dt(0, 1, df = 1), 0, )
              b[kk, z] ~ dnorm(0, sd = lambda[kk, z] * tau)
            }
          }
        } else {
          for (z in 1:(Z - 1)) {
            lambda[z] ~ T(dt(0, 1, df = 1), 0, )
            b[z] ~ dnorm(0, sd = lambda[z] * tau)
          }
        }

        if (Mlen > 1) {
          if (deg == 1) {
            for (z in 1:(Z - 1)) {
              for (zz in 1:Mlen) {
                lambda_m[zz, z] ~ T(dt(0, 1, df = 1), 0, )
                bm[zz, z] ~ dnorm(0, sd = lambda_m[zz, z] * tau)
              }
            }
          } else if (deg > 1) {
            for (z in 1:(Z - 1)) {
              for (zz in 1:Mlen) {
                # horseshoe prior
                lambda_m1[zz, z] ~ T(dt(0, 1, df = 1), 0, )
                bm1[zz, z] ~ dnorm(0, sd = lambda_m1[zz, z] * tau)
                lambda_m2[zz, z] ~ T(dt(0, 1, df = 1), 0, )
                bm2[zz, z] ~ dnorm(0, sd = lambda_m2[zz, z] * tau)
              }
            }
          }
        } else if (Mlen == 1) {
          if (deg == 1) {
            for (z in 1:(Z - 1)) {
              lambda_m[z] ~ T(dt(0, 1, df = 1), 0, )
              bm[z] ~ dnorm(0, sd = lambda_m[z] * tau)
            }
          } else if (deg > 1) {
            for (z in 1:(Z - 1)) {
              lambda_m1[z] ~ T(dt(0, 1, df = 1), 0, )
              bm1[z] ~ dnorm(0, sd = lambda_m1[z] * tau)
              lambda_m2[z] ~ T(dt(0, 1, df = 1), 0, )
              bm2[z] ~ dnorm(0, sd = lambda_m2[z] * tau)
            }
          }
        }
      } else if (ERRprior == "normal") {
        if (deg > 1) {
          for (z in 1:(Z - 1)) {
            for (k in 1:deg) {
              b[k, z] ~ dnorm(0, 0.001)
            }
          }
        } else {
          for (z in 1:(Z - 1)) {
            b[z] ~ dnorm(0, .001)
          }
        }

        if (Mlen > 1) {
          if (deg == 1) {
            for (z in 1:(Z - 1)) {
              for (k in 1:Mlen) {
                bm[k, z] ~ dnorm(0, 0.001)
              }
            }
          } else if (deg > 1) {
            for (z in 1:(Z - 1)) {
              for (k in 1:Mlen) {
                bm1[k, z] ~ dnorm(0, 0.001)
                bm2[k, z] ~ dnorm(0, 0.001)
              }
            }
          }
        } else if (Mlen == 1) {
          if (deg == 1) {
            for (z in 1:(Z - 1)) {
              bm[z] ~ dnorm(0, 0.001)
            }
          } else if (deg > 1) {
            for (z in 1:(Z - 1)) {
              bm1[z] ~ dnorm(0, 0.001)
              bm2[z] ~ dnorm(0, 0.001)
            }
          }
        }
      } else if (ERRprior == "doubleexponential") {
        if (deg > 1) {
          for (z in 1:(Z - 1)) {
            for (kk in 1:2) {
              lambda[kk, z] ~ T(dt(0, 1, df = 1), 0, )
              b[kk, z] ~ ddexp(0.0, lambda[kk, z])
            }
          }
        } else {
          for (z in 1:(Z - 1)) {
            lambda[z] ~ T(dt(0, 1, df = 1), 0, )
            b[z] ~ ddexp(0.0, lambda[z])
          }
        }

        if (Mlen > 1) {
          if (deg == 1) {
            for (z in 1:(Z - 1)) {
              for (zz in 1:Mlen) {
                lambda_m[zz, z] ~ T(dt(0, 1, df = 1), 0, )
                bm[zz, z] ~ ddexp(0.0, lambda_m[zz, z])
              }
            }
          } else if (deg > 1) {
            for (z in 1:(Z - 1)) {
              for (zz in 1:Mlen) {
                lambda_m1[zz, z] ~ T(dt(0, 1, df = 1), 0, )
                bm1[zz, z] ~ ddexp(0.0, lambda_m1[zz, z])
                lambda_m2[zz, z] ~ T(dt(0, 1, df = 1), 0, )
                bm2[zz, z] ~ ddexp(0.0, lambda_m2[zz, z])
              }
            }
          }
        } else if (Mlen == 1) {
          if (deg == 1) {
            for (z in 1:(Z - 1)) {
              lambda_m[z] ~ T(dt(0, 1, df = 1), 0, )
              bm[z] ~ ddexp(0.0, lambda_m[z])
            }
          } else if (deg > 1) {
            for (z in 1:(Z - 1)) {
              lambda_m1[z] ~ T(dt(0, 1, df = 1), 0, )
              bm1[z] ~ ddexp(0.0, lambda_m1[z])
              lambda_m2[z] ~ T(dt(0, 1, df = 1), 0, )
              bm2[z] ~ ddexp(0.0, lambda_m2[z])
            }
          }
        }
      }
    } else if (doseRRmod == "LINEXP") {
      if (ERRprior == "truncated_horseshoe") {
        tau ~ T(dt(0, 1, df = 1), 0, )
        for (z in 1:(Z - 1)) {
          # horseshoe prior
          lambda[1, z] ~ T(dt(0, 1, df = 1), 0, )
          b[1, z] ~ T(dnorm(0, sd = lambda[1, z] * tau), 0, )
          lambda[2, z] ~ T(dt(0, 1, df = 1), 0, )
          b[2, z] ~ dnorm(0, sd = lambda[2, z] * tau)
        }

        if (Mlen > 1) {
          for (z in 1:(Z - 1)) {
            for (zz in 1:Mlen) {
              # horseshoe prior
              lambda_m1[zz, z] ~ T(dt(0, 1, df = 1), 0, )
              bm1[zz, z] ~ T(dnorm(0, sd = lambda_m1[zz, z] * tau), 0, )
              lambda_m2[zz, z] ~ T(dt(0, 1, df = 1), 0, )
              bm2[zz, z] ~ dnorm(0, sd = lambda_m2[zz, z] * tau)
            }
          }
        } else if (Mlen == 1) {
          for (z in 1:(Z - 1)) {
            lambda_m1[z] ~ T(dt(0, 1, df = 1), 0, )
            bm1[z] ~ T(dnorm(0, sd = lambda_m1[z] * tau), 0, )
            lambda_m2[z] ~ T(dt(0, 1, df = 1), 0, )
            bm2[z] ~ dnorm(0, sd = lambda_m2[z] * tau)
          }
        }
      } else if (ERRprior == "truncated_normal") {
        for (z in 1:(Z - 1)) {
          b[1, z] ~ T(dnorm(0, 0.001), 0, )
          b[2, z] ~ dnorm(0, 0.001)
        }
        if (Mlen > 1) {
          for (z in 1:(Z - 1)) {
            for (k in 1:Mlen) {
              bm1[k, z] ~ T(dnorm(0, 0.001), 0, )
              bm2[k, z] ~ dnorm(0, 0.001)
            }
          }
        } else if (Mlen == 1) {
          for (z in 1:(Z - 1)) {
            bm1[z] ~ T(dnorm(0, 0.001), 0, )
            bm2[z] ~ dnorm(0, 0.001)
          }
        }
      } else if (ERRprior == "truncated_doubleexponential") {
        for (z in 1:(Z - 1)) {
          lambda[1, z] ~ T(dt(0, 1, df = 1), 0, )
          b[1, z] ~ T(ddexp(0.0, lambda[1, z]), 0, )
          lambda[2, z] ~ T(dt(0, 1, df = 1), 0, )
          b[2, z] ~ ddexp(0.0, lambda[2, z])
        }
        if (Mlen > 1) {
          for (z in 1:(Z - 1)) {
            for (zz in 1:Mlen) {
              lambda_m1[zz, z] ~ T(dt(0, 1, df = 1), 0, )
              bm1[zz, z] ~ T(ddexp(0.0, lambda_m1[zz, z]), 0, )
              lambda_m2[zz, z] ~ T(dt(0, 1, df = 1), 0, )
              bm2[zz, z] ~ ddexp(0.0, lambda_m2[zz, z])
            }
          }
        } else if (Mlen == 1) {
          for (z in 1:(Z - 1)) {
            lambda_m1[z] ~ T(dt(0, 1, df = 1), 0, )
            bm1[z] ~ T(ddexp(0.0, lambda_m1[z]), 0, )
            lambda_m2[z] ~ T(dt(0, 1, df = 1), 0, )
            bm2[z] ~ ddexp(0.0, lambda_m2[z])
          }
        }
      } else if (ERRprior == "horseshoe") {
        tau ~ T(dt(0, 1, df = 1), 0, )
        for (z in 1:(Z - 1)) {
          for (kk in 1:2) {
            # horseshoe prior
            lambda[kk, z] ~ T(dt(0, 1, df = 1), 0, )
            b[kk, z] ~ dnorm(0, sd = lambda[kk, z] * tau)
          }
        }

        if (Mlen > 1) {
          for (z in 1:(Z - 1)) {
            for (zz in 1:Mlen) {
              # horseshoe prior
              lambda_m1[zz, z] ~ T(dt(0, 1, df = 1), 0, )
              bm1[zz, z] ~ dnorm(0, sd = lambda_m1[zz, z] * tau)
              lambda_m2[zz, z] ~ T(dt(0, 1, df = 1), 0, )
              bm2[zz, z] ~ dnorm(0, sd = lambda_m2[zz, z] * tau)
            }
          }
        } else if (Mlen == 1) {
          for (z in 1:(Z - 1)) {
            lambda_m1[z] ~ T(dt(0, 1, df = 1), 0, )
            bm1[z] ~ dnorm(0, sd = lambda_m1[z] * tau)
            lambda_m2[z] ~ T(dt(0, 1, df = 1), 0, )
            bm2[z] ~ dnorm(0, sd = lambda_m2[z] * tau)
          }
        }
      } else if (ERRprior == "normal") {
        for (z in 1:(Z - 1)) {
          for (k in 1:2) {
            b[k, z] ~ dnorm(0, 0.001)
          }
        }

        if (Mlen > 1) {
          for (z in 1:(Z - 1)) {
            for (k in 1:Mlen) {
              bm1[k, z] ~ dnorm(0, 0.001)
              bm2[k, z] ~ dnorm(0, 0.001)
            }
          }
        } else if (Mlen == 1) {
          for (z in 1:(Z - 1)) {
            bm1[z] ~ dnorm(0, 0.001)
            bm2[z] ~ dnorm(0, 0.001)
          }
        }
      } else if (ERRprior == "doubleexponential") {
        for (z in 1:(Z - 1)) {
          for (kk in 1:2) {
            lambda[kk, z] ~ T(dt(0, 1, df = 1), 0, )
            b[kk, z] ~ ddexp(0.0, lambda[kk, z])
          }
        }

        if (Mlen > 1) {
          for (z in 1:(Z - 1)) {
            for (zz in 1:Mlen) {
              lambda_m1[zz, z] ~ T(dt(0, 1, df = 1), 0, )
              bm1[zz, z] ~ ddexp(0.0, lambda_m1[zz, z])
              lambda_m2[zz, z] ~ T(dt(0, 1, df = 1), 0, )
              bm2[zz, z] ~ ddexp(0.0, lambda_m2[zz, z])
            }
          }
        } else if (Mlen == 1) {
          for (z in 1:(Z - 1)) {
            lambda_m1[z] ~ T(dt(0, 1, df = 1), 0, )
            bm1[z] ~ ddexp(0.0, lambda_m1[z])
            lambda_m2[z] ~ T(dt(0, 1, df = 1), 0, )
            bm2[z] ~ ddexp(0.0, lambda_m2[z])
          }
        }
      }
    }
  }

  if (family == "gaussian") {
    sigma ~ T(dnorm(0, 0.001), 0, )
  }
  if (family == "prophaz") {
    for (j in 1:prophaz_numints) {
      h0[j] ~ dgamma(0.01, 0.01)
    }
  }

  # Likelihood

  if (family == "prophaz") {
    for (i in 1:N) {
      # Pieces of the cumulative baseline hazard function
      for (k in int.entry[i]:int.exit[i]) {
        start_t[i, k] <- max(prophaz_timepoints[k], entry[i])
        end_t[i, k] <- min(prophaz_timepoints[k + 1], exit[i])
        HH[i, k] <- max(end_t[i, k] - start_t[i, k], 0) * h0[k]
      }
      # Cumulative baseline hazard
      H[i] <- sum(HH[i, int.entry[i]:int.exit[i]])
    }
  }

  for (i in 1:N) {
    dose[i] <- inprod(dosemat[i, 1:K], vec[1:K])
    if (family != "multinomial") {
      if (doseRRmod != "LINEXP") {
        if (Mlen == 0) {
          if (deg == 2) {
            dosepart[i] <- b[1] * dose[i] + b[2] * dose[i] * dose[i]
          } else if (deg == 1) {
            dosepart[i] <- b * dose[i]
          }
        } else if (Mlen == 1) {
          if (deg == 2) {
            dosepart[i] <- (b[1] + bm1 * Mmat[i]) *
              dose[i] +
              (b[2] + bm2 * Mmat[i]) * dose[i] * dose[i]
          } else if (deg == 1) {
            dosepart[i] <- (b + bm * Mmat[i]) * dose[i]
          }
        } else {
          if (deg == 2) {
            dosepart[i] <- (b[1] + inprod(Mmat[i, ], bm1[1:Mlen])) *
              dose[i] +
              (b[2] + inprod(Mmat[i, ], bm2[1:Mlen])) * dose[i] * dose[i]
          } else if (deg == 1) {
            dosepart[i] <- (b + inprod(Mmat[i, ], bm[1:Mlen])) * dose[i]
          }
        }
      } else {
        # LINEXP
        if (Mlen == 0) {
          dosepart[i] <- b[1] * dose[i] * exp(b[2] * dose[i])
        } else if (Mlen == 1) {
          dosepart[i] <- (b[1] + bm1 * Mmat[i]) *
            dose[i] *
            exp((b[2] + bm2 * Mmat[i]) * dose[i])
        } else {
          dosepart[i] <- (b[1] + inprod(Mmat[i, ], bm1[1:Mlen])) *
            dose[i] *
            exp((b[2] + inprod(Mmat[i, ], bm2[1:Mlen])) * dose[i])
        }
      }

      if (Xlen > 1) {
        Xlinpred[i] <- a0 + inprod(Xmat[i, ], a[1:Xlen])
      } else if (Xlen == 1) {
        Xlinpred[i] <- a0 + Xmat[i] * a
      } else if (Xlen == 0) {
        Xlinpred[i] <- a0
      }
      if (family == "gaussian") {
        mu[i] <- dosepart[i] + Xlinpred[i]
        Y[i] ~ dnorm(mu[i], sd = sigma)
      } else if (family == "binomial") {
        if (doseRRmod %in% c("ERR", "LINEXP")) {
          mu[i] <- (1 + dosepart[i]) * exp(Xlinpred[i])
        } else if (doseRRmod == "EXP") {
          mu[i] <- exp(dosepart[i] + Xlinpred[i])
        }
        prob[i] <- mu[i] / (1 + mu[i])
        Y[i] ~ dbinom(prob[i], size = 1)
      } else if (family == "poisson") {
        if (doseRRmod %in% c("ERR", "LINEXP")) {
          mu[i] <- (1 + dosepart[i]) * exp(Xlinpred[i]) * P[i]
        } else if (doseRRmod == "EXP") {
          mu[i] <- exp(dosepart[i] + Xlinpred[i]) * P[i]
        }
        Y[i] ~ dpois(mu[i])
      } else if (family == "prophaz") {
        # Relative risk
        if (doseRRmod %in% c("ERR", "LINEXP")) {
          R[i] <- (1 + dosepart[i]) * exp(Xlinpred[i])
        } else if (doseRRmod == "EXP") {
          R[i] <- exp(dosepart[i] + Xlinpred[i])
        }

        # Log-hazard
        logHaz[i] <- log(h0[int.exit[i]] * R[i])
        # Log-survival
        logSurv[i] <- -H[i] * R[i]

        # Using the zeros trick
        phi[i] <- 100000 - delta[i] * logHaz[i] - logSurv[i]
        zeros[i] ~ dpois(phi[i])
      } else if (family == "clogit") {
        # Relative risk
        if (doseRRmod %in% c("ERR", "LINEXP")) {
          R[i] <- (1 + dosepart[i]) * exp(Xlinpred[i])
        } else if (doseRRmod == "EXP") {
          R[i] <- exp(dosepart[i] + Xlinpred[i])
        }
      }
    } else {
      # multinomial
      for (z in 1:(Z - 1)) {
        if (doseRRmod != "LINEXP") {
          if (Mlen == 0) {
            if (deg == 2) {
              dosepart[i, z] <- b[1, z] * dose[i] + b[2, z] * dose[i] * dose[i]
            } else if (deg == 1) {
              dosepart[i, z] <- b[z] * dose[i]
            }
          } else if (Mlen == 1) {
            if (deg == 2) {
              dosepart[i, z] <- (b[1, z] + bm1[z] * Mmat[i]) *
                dose[i] +
                (b[2, z] + bm2[z] * Mmat[i]) * dose[i] * dose[i]
            } else if (deg == 1) {
              dosepart[i, z] <- (b[z] + bm[z] * Mmat[i]) * dose[i]
            }
          } else {
            if (deg == 2) {
              dosepart[i, z] <- (b[1, z] + inprod(Mmat[i, ], bm1[1:Mlen, z])) *
                dose[i] +
                (b[2, z] + inprod(Mmat[i, ], bm2[1:Mlen, z])) *
                  dose[i] *
                  dose[i]
            } else if (deg == 1) {
              dosepart[i, z] <- (b[z] + inprod(Mmat[i, ], bm[1:Mlen, z])) *
                dose[i]
            }
          }
        } else {
          # LINEXP
          if (Mlen == 0) {
            dosepart[i, z] <- b[1, z] * dose[i] * exp(b[2, z] * dose[i])
          } else if (Mlen == 1) {
            dosepart[i, z] <- (b[1, z] + bm1[z] * Mmat[i]) *
              dose[i] *
              exp((b[2, z] + bm2[z] * Mmat[i]) * dose[i])
          } else {
            dosepart[i, z] <- (b[1, z] + inprod(Mmat[i, ], bm1[1:Mlen, z])) *
              dose[i] *
              exp((b[2, z] + inprod(Mmat[i, ], bm2[1:Mlen, z])) * dose[i])
          }
        }

        if (Xlen > 1) {
          Xlinpred[i, z] <- a0[z] + inprod(Xmat[i, ], a[1:Xlen, z])
        } else if (Xlen == 1) {
          Xlinpred[i, z] <- a0[z] + Xmat[i] * a[z]
        } else if (Xlen == 0) {
          Xlinpred[i, z] <- a0[z]
        }

        if (doseRRmod %in% c("ERR", "LINEXP")) {
          mu[i, z] <- (1 + dosepart[i, z]) * exp(Xlinpred[i, z])
        } else if (doseRRmod == "EXP") {
          mu[i, z] <- exp(dosepart[i, z] + Xlinpred[i, z])
        }
      }

      mu[i, Z] <- 1
      probs[i, 1:Z] <- mu[i, 1:Z] / sum(mu[i, 1:Z])
      Y[i, 1:Z] ~ dmulti(prob = probs[i, 1:Z], size = 1)
    }
  }

  if (family == "clogit") {
    for (setnr in 1:nsets) {
      #probs[setnr, 1:N] <- R[1:N]*setmat[setnr,1:N]/inprod(setmat[setnr, 1:N], R[1:N])
      #Ymat[setnr, 1:N] ~ dmulti(prob=probs[setnr, 1:N], size=1)
      probs[set_start[setnr]:set_end[setnr]] <- R[
        set_start[setnr]:set_end[setnr]
      ] /
        sum(R[set_start[setnr]:set_end[setnr]])
      Y[set_start[setnr]:set_end[setnr]] ~ dmulti(
        probs[set_start[setnr]:set_end[setnr]],
        size = 1
      )
    }
  }
})
