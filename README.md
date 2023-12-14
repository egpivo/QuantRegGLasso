## QuantRegGroups Package
[![R build status](https://github.com/egpivo/QuantRegGroups/workflows/R-CMD-check/badge.svg)](https://github.com/egpivo/QuantRegGroups/actions)
[![Coverage Status](https://img.shields.io/codecov/c/github/egpivo/QuantRegGroups/master.svg)](https://app.codecov.io/github/egpivo/QuantRegGroups?branch=master)

## Description
**QuantRegGroups** is an R package meticulously crafted to address adaptively weighted group Lasso procedures of quantile regression problems. It excels in simultaneous variable selection and structure identification for varying coefficient quantile regression models, as well as additive quantile regression models featuring ultra-high dimensional covariates.

## Installation
- **Install the current development version from GitHub**:
   ```r
   remotes::install_github("egpivo/QuantRegGroups")
   ```

**Please Note:**
- **Windows Users:** Ensure that you have [Rtools](https://cran.r-project.org/bin/windows/Rtools/) installed before proceeding with the installation.

- **Mac Users:** You need Xcode Command Line Tools and should install the library [`gfortran`](https://github.com/fxcoudert/gfortran-for-macOS/releases). Follow these steps in the terminal:
    ```bash
    brew update
    brew install gcc
    ```
    For a detailed solution, refer to this [link](https://thecoatlessprofessor.com/programming/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error/), or download and install the library [`gfortran`](https://github.com/fxcoudert/gfortran-for-macOS/releases) to resolve the "`ld: library not found for -lgfortran`" error.


### Author
- [Wen-Ting Wang](https://www.linkedin.com/in/wen-ting-wang-6083a17b)
- [Wei-Ying_Wu](https://projecteuclid.org/search?author=Wei-Ying_Wu)

 
### Maintainer
[Wen-Ting Wang](https://www.linkedin.com/in/wen-ting-wang-6083a17b)

### Reference
Toshio Honda, Ching-Kang Ing, Wei-Ying Wu (2019). [Adaptively weighted group Lasso for semiparametric quantile regression models](https://projecteuclid.org/journals/bernoulli/volume-25/issue-4B/Adaptively-weighted-group-Lasso-for-semiparametric-quantile-regression-models/10.3150/18-BEJ1091.full)").

## License
GPL-3

