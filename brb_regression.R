pooled_t_blb_multi <- function(
  formula, data, b = NULL, s = 20, r = 50,
  robust = TRUE, ridge = 1e-10, seed = NULL,
  contrasts.arg = NULL
) {
  # Returns: list(est, ci, pooled_q, se_global, settings)
  if (!is.null(seed)) set.seed(seed)
  stopifnot(is.data.frame(data), inherits(formula, "formula"))
  # --- terms/xlev from full data to fix factor levels & design columns ---
  terms_obj <- terms(formula, data = data)
  mf_full <- model.frame(terms_obj, data = data, na.action = na.pass)
  xlev <- .getXlevels(terms_obj, mf_full) # levels for factor terms
  N <- nrow(mf_full)
  if (is.null(b)) b <- max(1000L, as.integer(N^0.7))
  # Safe inverse helper: chol2inv → MASS::ginv → eigen pinv
  safe_inv <- function(A, ridge = 0, tol = 1e-12) {
    A2 <- A; diag(A2) <- diag(A2) + ridge
    out <- tryCatch(chol2inv(chol(A2)),
                    error = function(e) NULL)
    if (!is.null(out)) return(out)
    if (requireNamespace("MASS", quietly = TRUE)) {
      return(MASS::ginv(A2))
    }
    ev <- eigen(A2, symmetric = TRUE)
    nz <- ev$values > (max(ev$values, 0) * tol)
    if (!any(nz)) stop("Matrix is numerically singular.")
    V <- ev$vectors[, nz, drop = FALSE]
    Dinv <- diag(1 / ev$values[nz], nrow = sum(nz))
    V %*% Dinv %*% t(V)
  }
  build_Xy <- function(df) {
    mf <- model.frame(terms_obj, data = df, na.action = na.omit, xlev = xlev)
    y <- model.response(mf)
    X <- model.matrix(terms_obj, data = mf, contrasts.arg = contrasts.arg)
    list(X = X, y = y)
  }
  # Discover coeff names/dimension p from one subsample
  idx0 <- sample.int(N, min(b, N), replace = FALSE)
  xy0 <- build_Xy(data[idx0, , drop = FALSE])
  p <- ncol(xy0$X)
  coef_names <- colnames(xy0$X)
  # Storage
  sub_est <- matrix(NA_real_, nrow = s, ncol = p,
                    dimnames = list(NULL, coef_names))
  all_t <- matrix(NA_real_, nrow = s*r, ncol = p,
                  dimnames = list(NULL, coef_names))
  all_se <- matrix(NA_real_, nrow = s*r, ncol = p,
                   dimnames = list(NULL, coef_names))
  row_id <- 1L
  for (i in seq_len(s)) {
    # ----- outer subsample (without replacement) -----
    idx <- sample.int(N, min(b, N), replace = FALSE)
    xy <- build_Xy(data[idx, , drop = FALSE])
    X <- xy$X; y <- xy$y
    if (is.null(X) || nrow(X) == 0) next
    # guard: need residual df > 0
    if (nrow(X) <= ncol(X)) next
    # unweighted subsample fit
    fit <- lm.fit(X, y)
    beta_sub <- as.numeric(coef(fit))
    names(beta_sub) <- coef_names
    sub_est[i, ] <- beta_sub
    # ----- inner weighted bootstraps -----
    for (j in seq_len(r)) {
      # multinomial weights summing to N
      w <- as.numeric(rmultinom(1, size = N, prob = rep(1 / nrow(X), nrow(X))))
      if (!any(w)) next
      wfit <- lm.wfit(X, y, w)
      beta_star <- as.numeric(coef(wfit))
      # residuals & Bread
      u <- y - as.numeric(X %*% beta_star)
      XtWX <- crossprod(X * sqrt(w))
      Bread_inv <- safe_inv(XtWX, ridge = ridge)
      # variance (classical vs HC1)
      if (!robust) {
        rss <- sum(w * u^2)
        dof <- sum(w) - ncol(X)
        sigma2 <- rss / max(dof, 1)
        V <- Bread_inv * sigma2
      } else {
        # Meat = X' diag(w * u^2) X
        Xw <- X * sqrt(w * (u^2))
        Meat <- crossprod(Xw)
        n_eff <- sum(w) # equals N for BLB inner
        HC1 <- n_eff / max(n_eff - ncol(X), 1)
        V <- Bread_inv %*% Meat %*% Bread_inv * HC1
      }
      se_star <- sqrt(pmax(diag(V), 0))
      t_star <- (beta_star - beta_sub) / se_star
      all_se[row_id, ] <- se_star
      all_t[row_id, ] <- t_star
      row_id <- row_id + 1L
    }
  }
  # Drop unfilled rows
  if (row_id <= s*r) {
    keep <- seq_len(row_id - 1L)
    all_t <- all_t[keep, , drop = FALSE]
    all_se <- all_se[keep, , drop = FALSE]
  }
  # Bagged point estimate
  blb_est <- colMeans(sub_est, na.rm = TRUE)
  # Pooled percentile-t quantiles (2.5%, 97.5%) and global SE (median)
  q <- apply(all_t, 2, quantile, c(0.025, 0.975), na.rm = TRUE)
  se_global <- apply(all_se, 2, median, na.rm = TRUE)
  # Back-transform around bagged center
  ci <- rbind(blb_est - q[2, ] * se_global,
              blb_est - q[1, ] * se_global)
  rownames(ci) <- c("lower", "upper")
  colnames(ci) <- coef_names
  list(est = blb_est, ci = ci, pooled_q = q, se_global = se_global,
       settings = list(b = b, s = s, r = r, robust = robust, ridge = ridge))
}