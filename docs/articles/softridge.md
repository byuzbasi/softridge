# Soft-thresholded ridge regression with softridge

Version 0.3.0 implements the SVD, ridge path, soft-thresholding, exact
Gaussian risk, SURE surfaces, and cross-validation/refit scoring in
RcppArmadillo. The public R interfaces and statistical tuning rules are
retained. Native kernels do not draw random numbers; fold generation
remains in R. The relative singular value cutoff remains `1e-10`. OpenMP
is disabled in the package kernels.

Numerical equivalence to the archived R implementation is tested to a
tolerance of `1e-8`; floating-point results need not be bitwise
identical. Selection is also checked on the validation fixtures. A small
timing example does not establish production-scale speed or scientific
superiority.

``` r

library(softridge)
set.seed(1)
n <- 60; p <- 300
X <- matrix(rnorm(n * p), n, p)
beta <- c(3, -2, 1.5, rep(0, p - 3))
y <- drop(X %*% beta + rnorm(n))

## Cross-validated fit
cv <- cv_hr(X, y, nfolds = 5)
cv$selected
```

    ##   [1]   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18
    ##  [19]  19  20  21  22  23  24  25  26  27  28  29  30  31  32  33  34  35  36
    ##  [37]  37  38  39  40  41  42  43  44  45  46  47  48  49  50  51  52  53  54
    ##  [55]  55  56  57  58  59  60  61  62  63  64  65  66  67  68  69  70  71  72
    ##  [73]  73  74  75  76  77  78  79  80  81  82  83  84  85  86  87  88  89  90
    ##  [91]  91  92  93  94  95  96  97  98  99 100 101 102 103 104 105 106 107 108
    ## [109] 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126
    ## [127] 127 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143 144
    ## [145] 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159 160 161 162
    ## [163] 163 164 165 166 167 168 169 170 171 172 173 174 175 176 177 178 179 180
    ## [181] 181 182 183 184 185 186 187 188 189 190 191 192 193 194 195 196 197 198
    ## [199] 199 200 201 202 203 204 205 206 207 208 209 210 211 212 213 214 215 216
    ## [217] 217 218 219 220 221 222 223 224 225 226 227 228 229 230 231 232 233 234
    ## [235] 235 236 237 238 239 240 241 242 243 244 245 246 247 248 249 250 251 252
    ## [253] 253 254 255 256 257 258 259 260 261 262 263 264 265 266 267 268 269 270
    ## [271] 271 272 273 274 275 276 277 278 279 280 281 282 283 284 285 286 287 288
    ## [289] 289 290 291 292 293 294 295 296 297 298 299 300

``` r

## Exact risk of the estimator for given tuning values
Xs <- scale(X)
hr_risk(Xs, beta * attr(Xs, "scaled:scale"), sigma = 1, k = 10, gamma = 0.05)$risk
```

    ## [1] 0.6175626

## Training least-squares rescaling

`cv_hr(X, y, scale = TRUE)` estimates a single nonnegative scale from
each CV training set after thresholding. The final model repeats that
operation on all training observations. The scale is
`max(0, crossprod(Xs %*% theta_hat, yc) / sum((Xs %*% theta_hat)^2))`,
where `Xs` and `yc` are the training-centered and, when requested,
scaled design and centered response. A zero fitted vector receives scale
one. Scaling is applied before coefficients are mapped to the original
units; the intercept is recomputed. This option and `refit = TRUE` are
alternatives.

``` r

set.seed(9)
X_demo <- matrix(rnorm(24 * 36), 24)
y_demo <- drop(X_demo[, 1:2] %*% c(2, -1) + rnorm(24))
g_demo <- hr_grid(X_demo, y_demo, nk = 3, ng = 5)
f_demo <- cv_hr(X_demo, y_demo, nfolds = 3, grid = g_demo, scale = TRUE)
f_demo$scale
```

    ## [1] 5.527165

``` r

head(predict(hr(X_demo, y_demo, f_demo$k, f_demo$gamma, scale = TRUE), X_demo))
```

    ## [1] -1.3353871 -1.4168784 -0.3093913 -0.5326695  0.6387946 -2.0246942

`hr_risk(..., scale = c)` and `hr_sure(..., scale = c)` instead take a
**fixed** positive numeric scale. The finite-sample formulas and degrees
of freedom are valid for that fixed scale. They do not establish exact
risk or SURE for a training-estimated scale or for parameters selected
by CV/SURE.
