[![License](https://img.shields.io/badge/license-GPL(>=2)-C11B17.svg)](http://www.gnu.org/licenses/gpl-2.0.html)


# WienR
*R* package for the calculation of the partial derivatives of the PDF and CDF of the Wiener drift diffusion model


# Description
Calculate the partial derivative of the probability density (PDF) oder cumulative distribution function (CDF) of the Wiener drift diffusion with respect to the first-passage time t (only for PDF), the upper barrier a, the drift rate v, and the relative starting point w. In addition the PDF and CDF themselves are also provided. All calculations are done by using the logarithm, which makes it more stable.


# Installation
First `devtools` needs to be installed and then this package can be installed by using the `install_github` function with the full GitHub path "RaphaelHartmann/WienR"

```
install.packages("devtools")
library(devtools)
install_github("RaphaelHartmann/WienR")
```


# Usage
Since there is already a package that proviides the PDF and CDF for the Wiener drift diffusion model we decided not to follow the convention for PDFs (`dnorm`, `dunif`, etc.) and CDFs (`pnorm`, `punif`, etc.). Instead the PDF is called by `WienerPDF` and the CDF with `WienerCDF`. The derivative functions are named with a leading `dt`, `da`, `dv`, or `dw` indicating the partial derivative with respect to the first-passage time, the upper barrier, the drift rate, or the relative starting point, respectively. E.g. `daWienerPDF` is the function for the derivative of the PDF with respect to the upper barrier. The gradient functions are named with a leading `grad`.


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
