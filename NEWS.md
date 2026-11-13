# ModalNet 0.3.1

* The null model (intercept and dispersion only) is now always among the
  candidates compared by the EBIC. The automatic grid is built from the initial
  pseudo-response, so its first point could already contain active covariates
  (or a numerically negligible one), in which case the null model could not be
  selected. When that happens the null fit is computed directly and placed at
  the top of the path. Fits with a user-supplied `lambda` are unchanged.
* `penalty.factor` is rescaled to sum to `p`, as `glmnet` does internally.
  Before, with weights on another scale (adaptive Lasso) the line search
  compared a penalty different from the one minimized by the inner solver,
  rejected valid steps and the path reached `maxit`.
* At a grid point where the dispersion stays at the boundary of its range for
  ten consecutive iterations (exponential limit, mode near zero) the iteration
  now stops and the fit is flagged as not converged, instead of running to
  `maxit`.

# ModalNet 0.3.0

* First CRAN release.
* Dashun Liu added as author.
* Documentation of every exported function and method now describes its return
  value and includes runnable examples; references carry DOIs.
* `verbose = TRUE` reports progress through `message()`.
* `print()` labels `alpha = 0` fits as ridge.

# ModalNet 0.2.1

* Fixed the scale of the ridge part of the penalty for `alpha < 1`. glmnet
  standardizes the Gaussian response and therefore divides its ridge term by the
  standard deviation of the pseudo-response; the inner call now compensates, so
  the penalized objective is exactly
  `lambda * sum(alpha * |b| + (1 - alpha) * b^2 / 2)` on the standardized
  columns. Before the fix the line search rejected the glmnet step at some grid
  points and ridge paths reached `maxit`. Lasso fits are unchanged.

# ModalNet 0.2.0

* Package renamed from ModalLasso to ModalNet; main function `modal_net()`.
  `modal_lasso()` remains as an alias.
* Elastic net and ridge penalties through `alpha`.
