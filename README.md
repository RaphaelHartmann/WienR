[![License](https://img.shields.io/badge/license-GPL(>=2)-C11B17.svg)](http://www.gnu.org/licenses/gpl-2.0.html)


# WienR
*R* package for the calculation of the partial derivatives of the PDF and CDF of the Wiener drift diffusion model


# Description
Calculate the partial derivative of the probability density (PDF) oder cumulative distribution function (CDF) of the Wiener drift diffusion with respect to the first-passage time t (only for PDF), the upper barrier a, the drift rate v, and the relative starting point w. In addition the PDF and CDF themselves are also provided. All calculations are done by using the logarithm, which makes it more stable.


# Installation

## Requirements:
For the installation there are a few requirements. First of all, the *R* version must be 4.0 or above (https://cran.r-project.org/). Second, use `update.packages()` to update your packages for compatibility with *R* version 4.0 or above. The argument `ask=FALSE` can be used to update all packages without asking. Third, Windows users might have to install Rtools40 (https://cran.r-project.org/bin/windows/Rtools/).

## Installation from GitHub:
First `devtools` needs to be installed and then `WienR` can be installed by using the `install_github()` function with the full GitHub path "RaphaelHartmann/WienR":

  ```
  install.packages("devtools")
  library(devtools)
  install_github("RaphaelHartmann/WienR")
  ```


# Usage
Since there is already a package that proviides the PDF and CDF for the Wiener drift diffusion model we decided not to follow the convention for PDFs (`dnorm`, `dunif`, etc.) and CDFs (`pnorm`, `punif`, etc.). Instead the PDF is called by `WienerPDF` and the CDF with `WienerCDF`. The derivative functions are named with a leading `dt`, `da`, `dv`, or `dw` indicating the partial derivative with respect to the first-passage time, the upper barrier, the drift rate, or the relative starting point, respectively. E.g. `daWienerPDF` is the function for the derivative of the PDF with respect to the upper barrier. The gradient functions are named with a leading `grad`.

There are five main arguments, all vectorized:

`t`: the first-passage time

`response`: "upper" or "lower"

`a`: the upper barrier

`v`: the drift rate

`w`: the relative starting point

The length of all arguments must match, except the length is one. In case one argument has length one it will be replicated to match the length of the others.

And three optional arguments:

`precision`: the precision with which the calculation is done. The calculated value is guaranteed to be closer to the true value than this value.

`K`: number of components evaluated from the infinite sums for the formulae.

`n.threads`: number of threads for multithreading.

If only neither `precision` nor `K` is used, then a default precision of 12e-10 is used to calculate the number of components that guarantee the precision. If `precision` is used but not `K` the same happens but with the provided precision value. If both are provided the number of components that guarantees the precision is calculated and used, except if `K` is larger. If only `K` is provided then this is used as fixed number of components.

# Examples
```
t <- .75
response <- "upper"
a <- 1
v <- 4
w <- .5

WienerPDF(t, response, a, v, w)
WienerCDF(t, response, a, v, w)

dtWienerPDF(t, response, a, v, w)
daWienerPDF(t, response, a, v, w)
dvWienerPDF(t, response, a, v, w)
dwWienerPDF(t, response, a, v, w)

daWienerCDF(t, response, a, v, w)
dvWienerCDF(t, response, a, v, w)
dwWienerCDF(t, response, a, v, w)

gradWienerPDF(t, response, a, v, w)
gradWienerCDF(t, response, a, v, w)
```
