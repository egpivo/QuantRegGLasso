## QuantRegGLasso 1.0.0 (Release Date: 2024-01-17)
### Overview 
In this release, we have deployed the package to CRAN following its standard procedures. The main features of this release include:

- `qrglasso`: This function allows for model quantile regression by adaptively weighted group Lasso.
   - `predict`: Generate estimations.
   - `plot.qrglasso`: Investigate BIC performance.
   - `plot.qrglasso.predict`: Visualize estimations.

- `orthogonize_bspline`: Orthogonalize B-splines using the built-in function `splines::bs`.
---


## QuantRegGLasso 0.5.0 (Release Date: 2024-01-11)
### Overview 
- Added a `plot.qrglasso` function for displaying BIC w.r.t. hyperparameters via `qrglasso` object.

---

## QuantRegGLasso 0.4.0 (Release Date: 2024-01-09)
### Overview 
- Added a `plot.qrglasso.predict` function for displaying top k coefficient functions via `qrglasso.predict` object.

---

## QuantRegGLasso 0.3.0 (Release Date: 2024-01-08)
### Overview 
- Added a `predict` function for estimating top k coefficient functions via `qrglasso` object
- Added a helper function to meet the pre-conditions of `predict`

---

## QuantRegGLasso 0.2.0 (Release Date: 2024-01-06)
### Overview 
- Added one-dimensional example code in the main function `qrglasso`.
- Improved code quality and increased code coverage for the helper function `orthogonize_bspline`.

---
## QuantRegGLasso 0.1.0 (Release date: 2024-01-05)
#### Overview 
- Added a helper function to orthogonize B-splines

---
