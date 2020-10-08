[![License](https://img.shields.io/badge/license-GPL(>=2)-C11B17.svg)](http://www.gnu.org/licenses/gpl-2.0.html)


# WienR
Calculation of the partial derivatives of the PDF and CDF of the Wiener diffusion model


# Description
Calculate the partial derivative of the probability density (PDF) oder cumulative distribution function (CDF) of the Wiener drift diffusion with respect to the first-passage time t (only for PDF), the upper barrier a, the drift rate v, and the relative starting point w. In addition the PDF and CDF themselves are also provided. All calculations are done by using the logarithm, which makes it more stable.


# Installation
First `devtools` needs to be installed and then this package can be installed by using the `install_github` function with the full GitHub path "RaphaelHartmann/WienR"

```
install.packages("devtools")
library(devtools)
install_github("RaphaelHartmann/WienR")
```
