## QuantRegGLasso: Adaptively Weighted Group Lasso for Semiparametric Quantile Regression Models

[![License](https://eddelbuettel.github.io/badges/GPL2+.svg)](https://www.gnu.org/licenses/gpl-2.0.html)
[![R build status](https://github.com/egpivo/QuantRegGLasso/workflows/R-CMD-check/badge.svg)](https://github.com/egpivo/QuantRegGLasso/actions)
[![Code Coverage](https://img.shields.io/codecov/c/github/egpivo/QuantRegGLasso/master.svg)](https://app.codecov.io/github/egpivo/QuantRegGLasso?branch=master)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/QuantRegGLasso)](https://CRAN.R-project.org/package=QuantRegGLasso)
[![Downloads (monthly)](https://cranlogs.r-pkg.org/badges/QuantRegGLasso?color=brightgreen)](https://www.r-pkg.org/pkg/QuantRegGLasso)
[![Downloads (total)](https://cranlogs.r-pkg.org/badges/grand-total/QuantRegGLasso?color=brightgreen)](https://www.r-pkg.org/pkg/QuantRegGLasso)
[![BEJ](https://img.shields.io/badge/Bernoulli-10.3150%2FBEJ1091-brightgreen)](https://doi.org/10.3150/18-BEJ1091)


**QuantRegGLasso** is an R package designed for adaptively weighted group Lasso procedures in quantile regression. It excels in simultaneous variable selection and structure identification for varying coefficient quantile regression models and additive quantile regression models with ultra-high dimensional covariates.


## Installation
- Install the current development version from GitHub:
  ```r
  remotes::install_github("egpivo/QuantRegGLasso")
  ```

**Please Note:**

- **Windows Users:** Ensure that you have [Rtools](https://cran.r-project.org/bin/windows/Rtools/) installed before proceeding with the installation.

- **Mac Users:** You need Xcode Command Line Tools and should install the library [`gfortran`](https://github.com/fxcoudert/gfortran-for-macOS/releases). Follow these steps in the terminal:
  ```bash
  brew update
  brew install gcc
  ```

  For a detailed solution, refer to this [link](https://thecoatlessprofessor.com/programming/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error/), or download and install the library [`gfortran`](https://github.com/fxcoudert/gfortran-for-macOS/releases) to resolve the "`ld: library not found for -lgfortran`" error.


### Authors
- [Wen-Ting Wang](https://www.linkedin.com/in/wen-ting-wang-6083a17b) ([GitHub](https://github.com/egpivo))
- [Wei-Ying Wu](https://projecteuclid.org/search?author=Wei-Ying_Wu)
- [Toshio Honda](https://www1.econ.hit-u.ac.jp/honda/e-honda.html)
- [Ching-Kang Ing](https://www.researchgate.net/profile/Ching-Kang-Ing)

 
### Maintainer
[Wen-Ting Wang](https://www.linkedin.com/in/wen-ting-wang-6083a17b) ([GitHub](https://github.com/egpivo))

### Reference
Toshio Honda, Ching-Kang Ing, Wei-Ying Wu (2019). [Adaptively weighted group Lasso for semiparametric quantile regression models](https://projecteuclid.org/journals/bernoulli/volume-25/issue-4B/Adaptively-weighted-group-Lasso-for-semiparametric-quantile-regression-models/10.3150/18-BEJ1091.full).

This paper introduces the adaptively weighted group Lasso procedure and its application to semiparametric quantile regression models. The methodology is grounded in a strong sparsity condition, establishing selection consistency under certain weight conditions.

## License
GPL (>= 2)
