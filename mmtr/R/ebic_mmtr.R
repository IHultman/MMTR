#' This function computes the EBIC for an MMTR model.
#'
#' @param mmtr_mod MMTR model object produced by calling the mmtr() function.
#' @param y Numeric M-vector of responses across all N groups, where M = M1 + M2 + ... + MN.
#' @param X (D+1)-mode fixed effects tensor of dimension P1xP2x...xPDxM, i.e. there is a D-mode
#'   P1xP2x...xPD fixed effects covariate tensor associated with each of the M responses.
#' @param Z (D+1)-mode random effects tensor of dimension Q1xQ2x...xQDxM, i.e. there is a D-mode
#'   Q1xQ2x...xQD random effects covariate tensor associated with each of the M responses.
#' @param group_ids Integer or factor M-vector indicating which group each observation belongs to.
#' @param n_reps Maximum number of iterations to run the mixed model trace regression AECM algorithm.
#'   Defaults to 1000.
#' @param conv_thrsh Convergence threshold used for early stopping of the mixed model trace regression
#'   AECM algorithm based on the change in the log-likelihood from one iteration to the next.
#'   Defaults to 1e-3.
#' @param n_proc Number of processors to use for parallel computation. Defaults to one.
#' @export
ebic_mmtr = function(
  mmtr_mod,
  y, X, Z,
  group_ids,
  n_reps=1,
  conv_thrsh=1e-3,
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
  #           n_reps: Maximum number of iterations to run the mixed model trace regression
  #                   AECM algorithm.
  #
  #       conv_thrsh: Convergence threshold used for early stopping of the mixed model trace
  #                   regression AECM algorithm based on the change in the log-likelihood
  #                   from one iteration to the next.
  #
  #           n_proc: Number of processors to use for parallel computation.

  stopifnot((typeof(y) == "double") || (typeof(y) == "integer") );
  stopifnot(typeof(group_ids) == "integer");

  y = as.vector(y);
  group_ids = as.factor(as.vector(group_ids) );

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

  L_list = mmtr_mod$L_hat_list;

  stopifnot(typeof(L_list) == "list");
  stopifnot(length(L_list) == t_dim);

  for (kx in 1:t_dim) {
    stopifnot(is.element("matrix", class(L_list[[kx]]) ));
    stopifnot(nrow(L_list[[kx]]) == q_dims[kx]);
  }

  Lt_Zt_mats = tucker_product_matrix(Zt_mats, L_list);

  # Track convergence data
  neg_ll_vals = rep(NA, n_reps);

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

  for (rx in 1:n_reps) {
    B_results = ebic_mmtr_update_b_par(
      y,
      Xt_mats,
      Lt_Zt_mats,
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

    L_list = ebic_mmtr_update_L_par(
      y_tilde,
      Zt_mats,
      L_list, tau2,
      grps_ixs_table,
      proc_ixs_table);

    if (is.null(L_list) ) {
      return(NA);
    }

    Lt_Zt_mats = tucker_product_matrix(Zt_mats, L_list);

    neg_ll_vals[rx] = neg_ll_mmtr_par(
      y_tilde,
      Lt_Zt_mats,
      tau2,
      grps_ixs_table,
      proc_ixs_table);

    if ((rx > 1) && (abs(neg_ll_vals[rx] - neg_ll_vals[rx-1]) < conv_thrsh) ) {
      break;
    }
  }

  s_dims = unlist(lapply(L_list, ncol) );

  n_nz_params = n_nz_b + sum(q_dims * s_dims) + 1;
  n_tot_params = P_dim + sum(q_dims * s_dims) + 1;

  (2 * neg_ll_vals[rx]) +
  (log(n_dim) * n_nz_params) +
  (2 * lchoose(n_tot_params, n_nz_params) )
};


ebic_mmtr_update_b_par = function(
  y,
  Xt_mats,
  Lt_Zt_mats,
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
  #           Lt_Zt_mats: Matrix whose columns are the Kronecker product L2 (x) L1
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

  #start_time = Sys.time();

  X_w_y_breve = foreach::foreach(px=1:n_proc, .combine="rbind") %dopar% {
    proc_ax = proc_ixs_table[1,px];
    proc_bx = proc_ixs_table[2,px];

    proc_obs_ax = grps_ixs_table[1,proc_ax];
    proc_obs_bx = grps_ixs_table[2,proc_bx];

    Lt_Zt_mats_px = Lt_Zt_mats[,proc_obs_ax:proc_obs_bx];

    if (!is.element("matrix", class(Lt_Zt_mats) ) || (nrow(Lt_Zt_mats) == 1) ) {
      Lt_Zt_mats_px = matrix(Lt_Zt_mats_px, 1);
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

      Lt_Zt_mats_ix = Lt_Zt_mats_px[,ax:bx];

      if (nrow(Lt_Zt_mats_px) == 1) {
        Lt_Zt_mats_ix = matrix(Lt_Zt_mats_ix, 1);
      }

      gamma_sqrt_inv_ix = backsolve(
        chol(crossprod(Lt_Zt_mats_ix) + I_mix),
        I_mix,
        transpose=TRUE);

      Xt_breve_px[,ax:bx] = tcrossprod(Xt_mats_px[,ax:bx], gamma_sqrt_inv_ix);
      y_breve_px[ax:bx] = gamma_sqrt_inv_ix %*% y_px[ax:bx];
    }

    cbind(t(Xt_breve_px), y_breve_px)
  };

  #end_time = difftime(Sys.time(), start_time, units="secs")[[1]];

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


ebic_mmtr_update_L_par = function(
  y_tilde,
  Zt_mats,
  L_list, tau2,
  grps_ixs_table,
  proc_ixs_table)
{
  #   Description:
  #
  #     Arguments:
  #              y_tilde: The full vector of differences between the response vector y
  #                       and the fixed effects fitted values Xb, i.e. y_tilde = y - Xb.
  #
  #              Zt_mats: Matrix whose columns are the vectorized random effects covariate
  #                       tensors.
  #
  #               L_list: List of D matrices, each of which is the square root covariance
  #                       matrix to one of the D mode's random effects.
  #
  #                 tau2: Updated estimate of tau2.
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
  #    L_list: Updated estimate to the random effects square root covariance matrices.

  n_dim = ncol(grps_ixs_table);
  q_dims = unlist(lapply(L_list, nrow) );
  Q_dim = prod(q_dims);
  s_dims = unlist(lapply(L_list, ncol) );
  S_dim = prod(s_dims);
  t_dim = length(L_list);

  n_proc = ncol(proc_ixs_table);

  # Check if either mode's square root covariance is zeroed out.
  if (any(unlist(lapply(L_list, function(L_kx) all(L_kx == 0) )))) {
    L_list[[1]] = matrix(0, q_dims[1], 1);
    L_list[[2]] = matrix(0, q_dims[2], 1);

    return(L_list);
  }

  #start_time = Sys.time();

  H_hat_w_g1 = foreach::foreach(px=1:n_proc, .combine="+") %dopar% {
    H_hat_w_g1_px = matrix(0, q_dims[1] * s_dims[1], (q_dims[1] * s_dims[1]) + 1);

    Kq1q2_ixs = get_commutation_ixs(q_dims);
    Ks1s2_ixs = get_commutation_ixs(s_dims);
    Ks2q1_ixs = get_commutation_ixs(c(s_dims[2], q_dims[1]) );

    proc_ax = proc_ixs_table[1,px];
    proc_bx = proc_ixs_table[2,px];

    proc_obs_ax = grps_ixs_table[1,proc_ax];
    proc_obs_bx = grps_ixs_table[2,proc_bx];

    Z1t_mats_px = Zt_mats[,proc_obs_ax:proc_obs_bx];
    y_tilde_px = y_tilde[proc_obs_ax:proc_obs_bx];
    vec_L2t_Z2t_mats_px = matrix(
      crossprod(L_list[[2]], matrix(Z1t_mats_px[Kq1q2_ixs,], q_dims[2]) ),
      q_dims[1] * s_dims[2]);

    # This is equal to the Kronecker product of [L2t (x) L1t] times Z1t_mats_px.
    Lt_Z1t_mats_px = matrix(
      crossprod(L_list[[1]], matrix(vec_L2t_Z2t_mats_px[Ks2q1_ixs,], q_dims[1]) ),
      S_dim);

    for (ix in proc_ax:proc_bx) {
      ax = grps_ixs_table[1,ix] - proc_obs_ax + 1;
      bx = grps_ixs_table[2,ix] - proc_obs_ax + 1;
      m_ix = bx - ax + 1;

      Lt_Z1t_mats_ix = Lt_Z1t_mats_px[,ax:bx];

      if (S_dim == 1) {
        Lt_Z1t_mats_ix = matrix(Lt_Z1t_mats_ix, 1);
      }

      # Compute the conditional expectations mu_hat and sigma_hat.
      Z1_mats_LLt_Z1t_mats_ix = crossprod(Lt_Z1t_mats_ix);
      mu1_hat_ix = (
        Lt_Z1t_mats_ix %*%
        solve(Z1_mats_LLt_Z1t_mats_ix + diag(1, m_ix) ) %*%
        y_tilde_px[ax:bx]);

      sigma1_hat_ix = tau2 * solve(tcrossprod(Lt_Z1t_mats_ix) + diag(1, S_dim) );
      gamma_hat_ix = sigma1_hat_ix + tcrossprod(mu1_hat_ix);

      if (S_dim > 1) {
        gamma_hat_ix = gamma_hat_ix[Ks1s2_ixs,][,Ks1s2_ixs];
      }

      gamma_sqrt_hat_ix = t(chol(gamma_hat_ix) );

      Z1_mats_L2_ix = t(matrix(vec_L2t_Z2t_mats_px[,ax:bx], s_dims[2]) );
      Z1_mats_L2_gamma_sqrt = Z1_mats_L2_ix %*% matrix(gamma_sqrt_hat_ix, s_dims[2]);

      order_ixs = get_commutation_ixs(c(m_ix, s_dims[1]) );
      n_ord_ixs = length(order_ixs);
      order_ixs = (
        rep(order_ixs, S_dim) +
        (n_ord_ixs * rep(0:(S_dim-1), each=n_ord_ixs) ));

      H_hat1_ix = tcrossprod(
        matrix(
          matrix(Z1_mats_L2_gamma_sqrt, q_dims[1])[,order_ixs],
          q_dims[1] * s_dims[1]) );

      g_hat1_ix = as.vector(
        crossprod(
          matrix(vec_L2t_Z2t_mats_px[,ax:bx] %*% y_tilde_px[ax:bx], s_dims[2]),
          matrix(mu1_hat_ix[Ks1s2_ixs,], s_dims[2]) ));

      H_hat_w_g1_px = H_hat_w_g1_px + cbind(H_hat1_ix, g_hat1_ix);
    }

    H_hat_w_g1_px
  };

  #end_time = difftime(Sys.time(), start_time, units="secs")[[1]];
  #cat("H1_hat compute time: ", end_time, " s\n");
  #start_time = Sys.time();

  ix_g1 = (q_dims[1] * s_dims[1]) + 1;
  H_hat1_svd = svd(H_hat_w_g1[,-ix_g1]);
  diag_cutoff = 1e-12;
  rnk = sum(H_hat1_svd$d > diag_cutoff);

  if (rnk < 1) {
    return(NULL);
  }

  L_list[[1]] = matrix(
    tcrossprod(
      H_hat1_svd$u[,1:rnk] %*% diag(x=(1 / H_hat1_svd$d[1:rnk]), rnk),
      matrix(H_hat1_svd$v[,1:rnk], ncol=rnk) ) %*%
    H_hat_w_g1[,ix_g1],
    q_dims[1],
    s_dims[1]);

  #end_time = difftime(Sys.time(), start_time, units="secs")[[1]];
  #cat("L1_hat update time: ", end_time, " s\n");
  #start_time = Sys.time();

  H_hat_w_g2 = foreach::foreach(px=1:n_proc, .combine="+") %dopar% {
    H_hat_w_g2_px = matrix(0, q_dims[2] * s_dims[2], (q_dims[2] * s_dims[2]) + 1);

    Ks1q2_ixs = get_commutation_ixs(c(s_dims[1], q_dims[2]) );
    Ks2s1_ixs = get_commutation_ixs(rev(s_dims) );

    proc_ax = proc_ixs_table[1,px];
    proc_bx = proc_ixs_table[2,px];

    proc_obs_ax = grps_ixs_table[1,proc_ax];
    proc_obs_bx = grps_ixs_table[2,proc_bx];

    Z1t_mats_px = Zt_mats[,proc_obs_ax:proc_obs_bx];
    y_tilde_px = y_tilde[proc_obs_ax:proc_obs_bx];
    vec_L1t_Z1t_mats_px = matrix(
      crossprod(L_list[[1]], matrix(Z1t_mats_px, q_dims[1]) ),
      q_dims[2] * s_dims[1]);

    # This is equal to the Kronecker product of [L2t (x) L1t] times Z1t_mats_px.
    Lt_Z1t_mats_px = matrix(
      crossprod(L_list[[2]], matrix(vec_L1t_Z1t_mats_px[Ks1q2_ixs,], q_dims[2]) ),
      S_dim);

    if (S_dim > 1) {
      Lt_Z1t_mats_px = Lt_Z1t_mats_px[Ks2s1_ixs,];
    }

    for (ix in proc_ax:proc_bx) {
      ax = grps_ixs_table[1,ix] - proc_obs_ax + 1;
      bx = grps_ixs_table[2,ix] - proc_obs_ax + 1;
      m_ix = bx - ax + 1;

      Lt_Z1t_mats_ix = Lt_Z1t_mats_px[,ax:bx];

      if (S_dim == 1) {
        Lt_Z1t_mats_ix = matrix(Lt_Z1t_mats_ix, 1);
      }

      # Compute the conditional expectations mu_hat and sigma_hat.
      Z1_mats_LLt_Z1t_mats_ix = crossprod(Lt_Z1t_mats_ix);
      mu1_hat_ix = (
        Lt_Z1t_mats_ix %*%
        solve(Z1_mats_LLt_Z1t_mats_ix + diag(1, m_ix) ) %*%
        y_tilde_px[ax:bx]);

      sigma1_hat_ix = tau2 * solve(tcrossprod(Lt_Z1t_mats_ix) + diag(1, S_dim) );
      gamma1_hat_ix = sigma1_hat_ix + tcrossprod(mu1_hat_ix);
      gamma1_sqrt_hat_ix = t(chol(gamma1_hat_ix) );

      Z2_mats_L1_ix = t(matrix(vec_L1t_Z1t_mats_px[,ax:bx], s_dims[1]) );
      Z2_mats_L1_gamma1_sqrt = Z2_mats_L1_ix %*% matrix(gamma1_sqrt_hat_ix, s_dims[1]);

      order_ixs = get_commutation_ixs(c(m_ix, s_dims[2]) );
      n_ord_ixs = length(order_ixs);
      order_ixs = (
        rep(order_ixs, S_dim) +
        (n_ord_ixs * rep(0:(S_dim-1), each=n_ord_ixs) ));

      H_hat2_ix = tcrossprod(
        matrix(
          matrix(Z2_mats_L1_gamma1_sqrt, q_dims[2])[,order_ixs],
          q_dims[2] * s_dims[2]) );

      g_hat2_ix = as.vector(
        crossprod(
          matrix(vec_L1t_Z1t_mats_px[,ax:bx] %*% y_tilde_px[ax:bx], s_dims[1]),
          matrix(mu1_hat_ix, s_dims[1]) ));

      H_hat_w_g2_px = H_hat_w_g2_px + cbind(H_hat2_ix, g_hat2_ix);
    }

    H_hat_w_g2_px
  };

  #end_time = difftime(Sys.time(), start_time, units="secs")[[1]];
  #cat("H2_hat compute time: ", end_time, " s\n");
  #start_time = Sys.time();

  ix_g2 = (q_dims[2] * s_dims[2]) + 1;
  H_hat2_svd = svd(H_hat_w_g2[,-ix_g2]);
  diag_cutoff = 1e-12;
  rnk = sum(H_hat2_svd$d > diag_cutoff);

  if (rnk < 1) {
    return(NULL);
  }

  L_list[[2]] = matrix(
    tcrossprod(
      H_hat2_svd$u[,1:rnk] %*% diag(x=(1 / H_hat2_svd$d[1:rnk]), rnk),
      matrix(H_hat2_svd$v[,1:rnk], ncol=rnk) ) %*%
    H_hat_w_g2[,ix_g2],
    q_dims[2],
    s_dims[2]);

  #end_time = difftime(Sys.time(), start_time, units="secs")[[1]];
  #cat("L2_hat update time: ", end_time, " s\n");

  L_list
};
