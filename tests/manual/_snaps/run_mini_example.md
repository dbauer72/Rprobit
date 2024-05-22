# test #1: simple cross sectional probit with three choices

    Code
      round(est[[1]], 2)
    Output
                                           mse
      true   -1.02 -0.08 -0.23 -0.82 0.77 0.00
      SJ     -1.08 -0.09 -0.26 -0.40 1.08 0.06
      TVBS   -1.08 -0.09 -0.26 -0.40 1.08 0.06
      TVBSv2 -1.08 -0.09 -0.26 -0.40 1.08 0.06
      ME     -1.08 -0.09 -0.26 -0.37 1.10 0.06

# test #2: same as test 1, but probit estimation

    Code
      round(est[[1]], 2)
    Output
                                           mse
      true   -1.02 -0.08 -0.23 -0.82 0.77 0.00
      SJ     -1.08 -0.09 -0.26 -0.40 1.08 0.06
      TVBS   -1.08 -0.09 -0.26 -0.40 1.08 0.06
      TVBSv2 -1.08 -0.09 -0.26 -0.40 1.08 0.06
      ME     -1.08 -0.09 -0.26 -0.37 1.10 0.06

# test #3: same as test 1, but panel with three choices per person

    Code
      round(est[[1]], 2)
    Output
                                          mse
      true   -1.02 -0.08 -0.23 -0.82 0.77   0
      SJ     -1.01 -0.07 -0.20 -0.78 0.77   0
      TVBS   -1.01 -0.07 -0.20 -0.78 0.77   0
      TVBSv2 -1.01 -0.07 -0.20 -0.78 0.77   0
      ME     -1.01 -0.07 -0.20 -0.79 0.79   0

# test #4: same as test 3, probit estimation

    Code
      round(est[[1]], 2)
    Output
                                          mse
      true   -1.02 -0.08 -0.23 -0.82 0.77   0
      SJ     -1.01 -0.07 -0.20 -0.78 0.77   0
      TVBS   -1.01 -0.07 -0.20 -0.78 0.77   0
      TVBSv2 -1.01 -0.07 -0.20 -0.78 0.77   0
      ME     -1.01 -0.07 -0.20 -0.79 0.79   0

# test #5: same as test 3, Tp 1-3 per decision, macml estimation

    Code
      round(est[[1]], 2)
    Output
                                           mse
      true   -1.02 -0.08 -0.23 -0.82 0.77 0.00
      SJ     -1.00 -0.09 -0.20 -0.60 0.90 0.01
      TVBS   -1.00 -0.09 -0.20 -0.60 0.90 0.01
      TVBSv2 -1.00 -0.09 -0.20 -0.60 0.90 0.01
      ME     -1.00 -0.09 -0.20 -0.56 0.94 0.02

# test #6: same as test 5, probit estimation

    Code
      round(est[[1]], 2)
    Output
                                          mse
      true   -1.02 -0.08 -0.23 -0.82 0.77   0
      SJ     -1.01 -0.07 -0.19 -0.84 0.72   0
      TVBS   -1.01 -0.07 -0.19 -0.80 0.75   0
      TVBSv2 -1.01 -0.07 -0.19 -0.84 0.72   0
      ME     -1.01 -0.07 -0.18 -0.83 0.75   0

# test #7: same as test 5, 2 random coefficients, macml estimation

    Code
      round(est[[1]], 2)
    Output
                                                         mse
      true -1.02 -0.08 -0.23 -0.82 0.77 -0.17 0.97 1.72 0.00
      SJ   -1.02 -0.06 -0.22 -0.63 0.83 -0.26 1.30 1.64 0.02

# test #8: same as test 7, probit  estimation

    Code
      round(est[[1]], 2)
    Output
                                                          mse
      true -1.02 -0.08 -0.23 -0.82 0.77 -0.17  0.97 1.72 0.00
      SJ   -1.07 -0.09 -0.20 -0.92 0.75 -0.18 -0.17 2.12 0.18

# test #9: same as test 7, only 1 random coefficient

    Code
      round(est[[1]], 2)
    Output
                                               mse
      true -1.02 -0.08 -0.23 -0.82 0.77 -0.17 0.00
      SJ   -0.96 -0.04  0.00 -0.84 0.64 -0.39 0.02

# test #10: same as test 9, probit estimation

    Code
      round(est[[1]], 2)
    Output
                                               mse
      true -1.02 -0.08 -0.23 -0.82 0.77 -0.17 0.00
      SJ   -0.96 -0.04  0.00 -0.84 0.64 -0.39 0.02

# test #11: V1 is alternative specific regressor, macml estimation

    Code
      round(est[[1]], 2)
    Output
                                                    mse
      true -1.02 -0.08 -0.23 -0.82 0.77 -0.17 0.97 0.00
      SJ   -1.07 -0.08 -0.21 -0.99 0.85 -0.14 1.00 0.01

# test #12: same as test 11, probit estimation

    Code
      round(est[[1]], 2)
    Output
                                                    mse
      true -1.02 -0.08 -0.23 -0.82 0.77 -0.17 0.97 0.00
      SJ   -1.07 -0.08 -0.21 -0.99 0.85 -0.14 1.00 0.01

# test #13: five categories, two regressors, one random coefficient

    Code
      round(est[[1]], 2)
    Output
                                              mse
      true -1.02 -0.08 -0.23 -0.82 0.77 -0.17   0
      SJ   -1.01 -0.06 -0.25 -0.77 0.76 -0.10   0

# test #15: ordered probit, estimated using the CML approach

    Code
      round(est[[1]], 2)
    Output
                                 mse
      true -1.02 1.47 -2.00 0.55   0
      SJ   -1.03 1.42 -2.08 0.59   0

# test #16: same as test 15, probit estimation

    Code
      round(est[[1]], 2)
    Output
                             mse
      true -1.02 1.47 -2.00 0.00
      SJ   -1.09 1.42 -2.11 0.01

# test #17: comparison to mlogit via Train dataset works

    Code
      round(comparison, 2)
    Output
                  price  time change comfort
      logit           0  0.03   0.33    0.95
      probit          0 -0.02  -0.22   -0.66
      logit_norm      1 19.32 219.85  637.12
      probit_norm     1 19.55 223.22  655.54

