[![License](https://img.shields.io/badge/license-GPL(>=2)-C11B17.svg)](http://www.gnu.org/licenses/gpl-2.0.html)


# WienR
*R* package for the calculation of the partial derivatives of the first-passage time PDF and CDF of the Wiener diffusion model with 3 to 7 parameters


# Description
Calculate the partial derivative of the first-passage time probability density (PDF) or cumulative distribution function (CDF) of the Wiener diffusion model with respect to the first-passage time `t` (only for PDF), the upper barrier `a`, the drift rate `v`, the relative starting point `w`, and in addition the non-decision time `t0`, the inter-trial variability of the drift rate `sv`, the inter-trial variability of the relative starting point `sw`, the inter-trial variability of the non-decision time `st0`. In addition the PDF and CDF themselves are also provided. Almost all calculations are done by using the logarithm, which makes it more stable.


# Installation

## Requirements:
For the installation there are a few requirements. First of all, the *R* version must be 4.0 or above (https://cran.r-project.org/). Second, use `update.packages()` to update your packages for compatibility with *R* version 4.0 or above. The argument `ask=FALSE` can be used to update all packages without asking each time. Third, Windows users might have to install Rtools40 (https://cran.r-project.org/bin/windows/Rtools/).

## Installation from GitHub:
First `devtools` needs to be installed and then `WienR` can be installed by using the `install_github()` function with the full GitHub path "RaphaelHartmann/WienR":

  ```
  install.packages("devtools")
  library(devtools)
  install_github("RaphaelHartmann/WienR")
  ```


# Usage
Since there is already a package that proviides the PDF and CDF for the Wiener diffusion model (with 4 parameters) we decided not to follow the convention for PDFs (`dnorm`, `dunif`, etc.) and CDFs (`pnorm`, `punif`, etc.). Instead the PDF is called by `WienerPDF` and the CDF with `WienerCDF`. The derivative functions are named with a leading `dt`, `da`, `dv`, `dw`, `dt0`, `dsv`, `dsw`, or `dst0` indicating the partial derivative with respect to the first-passage time, the upper barrier, the drift rate, the relative starting point, the non-decision time, the inter-trial variability of the drift rate, the inter-trial variability of the relative starting point, the inter-trial variability of the non-decision time, respectively. E.g. `daWienerPDF` is the function for the derivative of the PDF with respect to the upper barrier. The gradient functions are named with a leading `grad`.

There are nine main arguments, all vectorized:

`t`: the first-passage time

`response`: "upper" or "lower"

`a`: the upper barrier

`v`: the drift rate

`w`: the relative starting point

`t0`: the non-decision time

`sv`: the inter-trial variability of the drift rate

`sw`: the inter-trial variability of the relative starting point

`st0`: the inter-trial variability of the non-decision time


The length of all arguments must match, except if the length of one is used. In case one argument has length one it will be replicated to match the length of the others.

And three optional arguments:

`precision`: the precision with which the calculation is done. The calculated value is guaranteed to be closer to the true value than this value.

`K`: number of components evaluated from the infinite sums for the formulae.

`n.threads`: number of threads for multithreading.

If neither `precision` nor `K` is used, then a default precision of 1e-12 is used to calculate the number of components that guarantee the precision. If `precision` is used but not `K` the same happens but with the provided precision value. If both are provided the number of components that guarantees the precision is calculated and used, except if `K` is larger. If only `K` is provided then this is used as fixed number of components.

## Note
With the use of `sv`, `sw`, and/or `st0` numerical integration is needed. This is done by using the *C* library `cubature` by Steven G. Johnson (https://github.com/stevengj/cubature). This numerical integration allows for specifying the absolute and relative error tolerance. In the package we use the absolute error. By specifying `precision = 1e-10` the absolute error for the core function (PDF, CDF, or derivative function of the 3-parameter model) is set to 1e-11 and the absolute error for the numerical integral to 9e-11. Since the absolute error tolerances are additive, they sum to 1e-10. Using numerical integration always takes some time, especially with high dimensional integrals and a low absolute error. Therefore, the calculations of the PDF, CDF as well as their derivatives takes longer when using the additional parameters `sv`, `sw`, and/or `st0`.

# Examples
```
t <- .75
response <- "upper"
a <- 1
v <- .7
w <- t0 <- .5
sv <- sw <- st0 <- .1

WienerPDF(t, response, a, v, w)
WienerPDF(t, response, a, v, w, t0)
WienerPDF(t, response, a, v, w, t0, sv)
WienerPDF(t, response, a, v, w, t0, sv, sw)
WienerPDF(t, response, a, v, w, t0, sv, sw, st0)

WienerCDF(t, response, a, v, w, t0, st0)

dtWienerPDF(t, response, a, v, w)
daWienerPDF(t, response, a, v, w, t0)
dvWienerPDF(t, response, a, v, w, t0, sv)
dwWienerPDF(t, response, a, v, w, t0, sv, sw)
dt0WienerPDF(t, response, a, v, w, t0, sv, sw, st0)
dswWienerPDF(t, response, a, v, w, t0, sv, sw, st0)

dst0WienerCDF(t, response, a, v, w, t0, sv, sw, st0)

gradWienerPDF(t, response, a, v, w, t0, sv, sw, st0)
gradWienerCDF(t, response, a, v, w, t0, sv, sw, st0)
```
