#' This function computes the EBIC for an MMTR model.
#'
#' @param mmtr_mod MMTR model object produced by calling the mmtr() function.
#' @param y Numeric M-vector of responses across all N groups, where M = M1 + M2 + ... + MN.
#' @param X (D+1)-mode fixed effects tensor of dimension P1xP2x...xPDxM, i.e. there is a D-mode
#'   P1xP2x...xPD fixed effects covariate tensor associated with each of the M responses.
#' @param Z (D+1)-mode random effects tensor of dimension Q1xQ2x...xQDxM, i.e. there is a D-mode
#'   Q1xQ2x...xQD random effects covariate tensor associated with each of the M responses.
#' @param group_ids Integer or factor M-vector indicating which group each observation belongs to.
#' @param n_proc Number of processors to use for parallel computation. Defaults to one.
#' @export
bic_quasi_mmtr = function(
  mmtr_mod,
  y, X, Z,
  group_ids,
  n_proc=1)
{
  #   Description:
  #
  #     Arguments:
  #          mmtr_mod: MMTR model object produced by calling the mmtr() function.
  #
  #                 y: Numeric M-vector of responses across all N groups, where
  #                    M = M1 + M2 + ... + MN.
  #
  #                 X: (D+1)-mode fixed effects tensor of dimension P1xP2x...xPDxM, i.e.
  #                    there is a D-mode P1xP2x...xPD fixed effects covariate tensor
  #                    associated with each of the M responses.
  #
  #                 Z: (D+1)-mode random effects tensor of dimension Q1xQ2x...xQDxM, i.e.
  #                    there is a D-mode Q1xQ2x...xQD random effects covariate tensor
  #                    associated with each of the M responses.
  #
  #        group_ids: Integer or factor M-vector indicating which group each observation
  #                   belongs to.
  #
  #           n_proc: Number of processors to use for parallel computation.

  stopifnot((typeof(y) == "double") || (typeof(y) == "integer") );
  stopifnot(typeof(group_ids) == "integer");

  y = as.vector(y);
  group_ids = as.vector(group_ids);

  if (typeof(y) == "integer") {
    y = as.numeric(y);
  }

  stopifnot(is.element("array", class(X) ));
  stopifnot(is.element("array", class(Z) ));

  m_dim = length(y);
  n_dim = length(unique(group_ids) );
  t_dim = length(dim(X) ) - 1;

  stopifnot(t_dim == 2);
  stopifnot(length(group_ids) == m_dim);
  stopifnot(length(dim(Z) ) == (t_dim + 1) );
  stopifnot(dim(X)[t_dim+1] == m_dim);
  stopifnot(dim(Z)[t_dim+1] == m_dim);

  p_dims = dim(X)[1:t_dim];
  P_dim = prod(p_dims);
  q_dims = dim(Z)[1:t_dim];
  Q_dim = prod(q_dims);

  # Sort observations here so that group IDs are clustered.
  grps_order_ixs = order(group_ids);
  group_ids = group_ids[grps_order_ixs];
  y = y[grps_order_ixs];
  X = X[,,grps_order_ixs];
  Z = Z[,,grps_order_ixs];

  # Index ranges for each group.
  grps_info_table = as.data.frame(table(group_ids) );
  colnames(grps_info_table) = c("grp_id", "m_gx");
  grps_info_table$ax = cumsum(c(1, grps_info_table$m_gx[-n_dim]) );
  grps_info_table$bx = cumsum(grps_info_table$m_gx);

  # Check to make sure group ordering matches.
  for (gx in 1:n_dim) {
    grp_gx_ixs = grps_info_table$ax[gx]:grps_info_table$bx[gx];
    stopifnot(all(group_ids[grp_gx_ixs] == grps_info_table$grp_id[gx]) );
  }

  grps_ixs_table = t(as.matrix(grps_info_table[,c("ax", "bx")]) );
  rownames(grps_ixs_table) = NULL;

  Xt_mats = matrix(X, P_dim, m_dim);
  Zt_mats = matrix(Z, Q_dim, m_dim);

  stopifnot(is.element("matrix", class(mmtr_mod$B_hat) ));
  stopifnot(all(dim(mmtr_mod$B_hat) == p_dims) );

  ixs_b_nz = abs(as.vector(mmtr_mod$B_hat) ) > 0;
  n_nz_b = sum(ixs_b_nz);

  if (n_nz_b > 0) {
    Xt_mats = matrix(Xt_mats[ixs_b_nz,], n_nz_b);
  } else {
    Xt_mats = matrix(0, 1, m_dim);
  }

  if (n_dim < n_proc) {
    n_proc = n_dim;
  }

  grps_per_proc = floor(n_dim / n_proc);
  n_grps_leftover = n_dim %% n_proc;
  proc_sizes = rep(
    c(grps_per_proc + 1, grps_per_proc),
    c(n_grps_leftover, n_proc - n_grps_leftover) );

  proc_ixs_table = t(cbind(
    cumsum(c(1, proc_sizes[-n_proc]) ),
    cumsum(proc_sizes) ));

  doMC::registerDoMC(n_proc);

  B_results = bic_quasi_mmtr_update_b_par(
    y,
    Xt_mats,
    Zt_mats,
    grps_ixs_table,
    proc_ixs_table);

  b_vec = B_results$b_vec;
  tau2 = B_results$tau2;

  tau2_cutoff = 1e-8;
  if (tau2 < tau2_cutoff) {
    # Variance is fully explained by mean structure. EBIC cannot be computed.
    return(NA);
  }

  y_tilde = as.vector(y - crossprod(Xt_mats, b_vec) );

  quasi_neg_ll = neg_ll_mmtr_par(
    y_tilde,
    Zt_mats,
    tau2,
    grps_ixs_table,
    proc_ixs_table);

  s_dims = unlist(lapply(mmtr_mod$L_hat_list, ncol) );

  n_effective_pars = n_nz_b + sum(q_dims * s_dims) - prod(s_dims)^2;

  (2 * quasi_neg_ll) +
  (log(m_dim) * n_effective_pars)
};


bic_quasi_mmtr_update_b_par = function(
  y,
  Xt_mats,
  Zt_mats,
  grps_ixs_table,
  proc_ixs_table)
{
  #   Description:
  #
  #     Arguments:
  #                    y: Numeric M-vector of responses across all N groups, where
  #                       M = M1 + M2 + ... + MN.
  #
  #              Xt_mats: Matrix whose columns are the vectorized fixed effects covariate
  #                       tensors.
  #
  #           Zt_mats: Matrix whose columns are the Kronecker product L2 (x) L1
  #                       transposed times the vectorizations of each observation's Z_ij
  #                       matrix, i.e. the (i, j)th column is [L2^T (x) L1^T] vec(Z_ij).
  #
  #       grps_ixs_table: A table consisting of two rows, each column of which indicates
  #                       the first and last indices of the observations corresponding
  #                       to a particular group, i.e. the first value of the cth column
  #                       indicates the index of cth group's first observation, and the
  #                       second value of the cth column indicates the index of the cth
  #                       group's last observation.
  #
  #       proc_ixs_table: A table consisting of two rows, each column of which indicates
  #                       the first and last indices of the groups whose data will be
  #                       allocated to a particular processor for parallel computation.
  #
  #  Return value(s):
  #    b_vec: Updated estimate of the vectorized fixed effects tensor.
  #
  #     tau2: Updated estimate of tau2.

  n_dim = ncol(grps_ixs_table);
  m_dim = length(y);

  n_proc = ncol(proc_ixs_table);

  X_w_y_breve = foreach::foreach(px=1:n_proc, .combine="rbind") %dopar% {
    proc_ax = proc_ixs_table[1,px];
    proc_bx = proc_ixs_table[2,px];

    proc_obs_ax = grps_ixs_table[1,proc_ax];
    proc_obs_bx = grps_ixs_table[2,proc_bx];

    Zt_mats_px = Zt_mats[,proc_obs_ax:proc_obs_bx];

    if (!is.element("matrix", class(Zt_mats) ) || (nrow(Zt_mats) == 1) ) {
      Zt_mats_px = matrix(Zt_mats_px, 1);
    }

    Xt_mats_px = Xt_mats[,proc_obs_ax:proc_obs_bx];

    if (!is.element("matrix", class(Xt_mats) ) || (nrow(Xt_mats) == 1) ) {
      Xt_mats_px = matrix(Xt_mats_px, 1);
    }

    y_px = y[proc_obs_ax:proc_obs_bx];

    Xt_breve_px = array(NA, dim=dim(Xt_mats_px) );
    y_breve_px = rep(NA, length(y_px) );

    for (ix in proc_ax:proc_bx) {
      ax = grps_ixs_table[1,ix] - proc_obs_ax + 1;
      bx = grps_ixs_table[2,ix] - proc_obs_ax + 1;

      I_mix = diag(1, bx - ax + 1);

      Zt_mats_ix = Zt_mats_px[,ax:bx];

      if (nrow(Zt_mats_px) == 1) {
        Zt_mats_ix = matrix(Zt_mats_ix, 1);
      }

      gamma_sqrt_inv_ix = backsolve(
        chol(crossprod(Zt_mats_ix) + I_mix),
        I_mix,
        transpose=TRUE);

      Xt_breve_px[,ax:bx] = tcrossprod(Xt_mats_px[,ax:bx], gamma_sqrt_inv_ix);
      y_breve_px[ax:bx] = gamma_sqrt_inv_ix %*% y_px[ax:bx];
    }

    cbind(t(Xt_breve_px), y_breve_px)
  };

  ix_y_breve = ncol(X_w_y_breve);
  X_breve_svd = svd(X_w_y_breve[,-ix_y_breve]);
  diag_cutoff = 1e-12;
  rnk = sum(X_breve_svd$d > diag_cutoff);

  if (rnk < 1) {
    b_vec = matrix(0, nrow(Xt_mats) );
  } else {
    b_vec = (
      tcrossprod(
        X_breve_svd$v[,1:rnk] %*% diag(x=(1 / X_breve_svd$d[1:rnk]), rnk),
        matrix(X_breve_svd$u[,1:rnk], ncol=rnk) ) %*%
      X_w_y_breve[,ix_y_breve]);
  }

  tau2 = drop(
    crossprod(
      X_w_y_breve[,ix_y_breve] -
      (X_w_y_breve[,-ix_y_breve] %*% b_vec) )) / m_dim;

  list(b_vec=b_vec, tau2=tau2)
};
