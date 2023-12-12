## grpQuanReg Package

[![Travis-CI Build Status](https://travis-ci.org/egpivo/grpQuanReg.svg?branch=master)](https://travis-ci.org/egpivo/grpQuanReg)

## Description
**grpQuanReg** is an R package designed for meticulously crafted to address adaptively weighted group Lasso procedures. It excels in simultaneous variable selection and structure identification for varying coefficient quantile regression models, as well as additive quantile regression models featuring ultra-high dimensional covariates.

## Installation
There are two ways to install the package:
- **Install the current development version from GitHub**:
   ```r
   remotes::install_github("egpivo/grpQuanReg")
   ```

For compiling C++ code with the required [`RcppArmadillo`](https://CRAN.R-project.org/package=RcppArmadillo) and [`RcppParallel`](https://CRAN.R-project.org/package=RcppParallel) packages, follow these instructions:

* Windows users: Install [Rtools](https://CRAN.R-project.org/bin/windows/Rtools/)
* Mac users: Install Xcode Command Line Tools, and install the `gfortran` library. You can achieve this by running the following commands in the terminal:
  ```bash
  brew update
  brew install gcc
  ```

For a detailed solution, refer to [this link](https://thecoatlessprofessor.com/programming/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error/), or download and install the library [`gfortran`](https://github.com/fxcoudert/gfortran-for-macOS/releases) to resolve the error `ld: library not found for -lgfortran`.


### Author
- [Wen-Ting Wang](https://www.linkedin.com/in/wen-ting-wang-6083a17b)

 
### Maintainer
[Wen-Ting Wang](https://www.linkedin.com/in/wen-ting-wang-6083a17b)

### Reference
Toshio Honda, Ching-Kang Ing, Wei-Ying Wu (2019). [Adaptively weighted group Lasso for semiparametric quantile regression models](https://projecteuclid.org/journals/bernoulli/volume-25/issue-4B/Adaptively-weighted-group-Lasso-for-semiparametric-quantile-regression-models/10.3150/18-BEJ1091.full)").

## License
GPL-3

