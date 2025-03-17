get_commutation_ixs = function(dims){
  #   Description: Let X be a matrix with r rows and c columns. To get the order of indices
  #                that rearrange vec(X) to vec(X^T), form a matrix of integers 1:(r*c)
  #                with r rows and c columns, transpose the matrix and vectorize it.
  #                The resulting order will rearrange vec(X) to vec(X^T).
  #
  #     Arguments:
  #       dims: Vector whose first value is the number of rows and whose second value
  #             is the number of columns of a matrix X whose vectorization needs to be
  #             reordered to the vectorization of X^T.
  #
  #  Return value:

  as.vector(t(matrix(1:prod(dims), dims[1], dims[2]) ))
};


mu_hat_c_given_y = function(y_tilde, Lt_Zt_mats, grps_ixs_table) {
  #   Description: This function computes for each group the conditional expectation
  #                of that group's random effects vector c given y.
  #
  #     Arguments:
  #              y_tilde: The full vector of differences between the response vector y
  #                       and the fixed effects fitted values Xb, i.e. y_tilde = y - Xb.
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
  #  Return value: This function returns a matrix whose cth column is the conditional
  #                expectation of the cth group's random effects vector c_i given y_i.

  S_dim = nrow(Lt_Zt_mats);
  n_dim = ncol(grps_ixs_table);

  mu_hat_vecs = matrix(NA, S_dim, n_dim);

  for (gx in 1:n_dim) {
    grp_ax = grps_ixs_table[1,gx];
    grp_bx = grps_ixs_table[2,gx];
    m_gx = grp_bx - grp_ax + 1;

    Lt_Zt_mats_gx = Lt_Zt_mats[,grp_ax:grp_bx];

    if (S_dim == 1) {
      Lt_Zt_mats_gx = matrix(Lt_Zt_mats_gx, 1);
    }

    # Compute the conditional expectations mu_hat.
    Z_mats_LLt_Zt_mats_gx = crossprod(Lt_Zt_mats_gx);
    mu_hat_vecs[,gx] = as.matrix(
      Lt_Zt_mats_gx %*%
      solve(Z_mats_LLt_Zt_mats_gx + diag(1, m_gx) ) %*%
      y_tilde[grp_ax:grp_bx]);
  }

  mu_hat_vecs
};


neg_ll_mmtr_par = function(
  y_tilde,
  Lt_Zt_mats,
  tau2,
  grps_ixs_table,
  proc_ixs_table)
{
  #   Description:
  #
  #     Arguments:
  #              y_tilde: The full vector of differences between the response vector y
  #                       and the fixed effects fitted values Xb, i.e. y_tilde = y - Xb.
  #
  #           Lt_Zt_mats: Matrix whose columns are the Kronecker product L2 (x) L1
  #                       transposed times the vectorizations of each observation's Z_ij
  #                       matrix, i.e. the (i, j)th column is [L2^T (x) L1^T] vec(Z_ij).
  #
  #                 tau2: Updated estimate to tau2.
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
  #  Return value:

  m_dim = length(y_tilde);
  n_dim = ncol(grps_ixs_table);

  n_proc = ncol(proc_ixs_table);

  neg_ll = foreach::foreach(px=1:n_proc, .combine="+") %dopar% {
    proc_ax = proc_ixs_table[1,px];
    proc_bx = proc_ixs_table[2,px];

    proc_obs_ax = grps_ixs_table[1,proc_ax];
    proc_obs_bx = grps_ixs_table[2,proc_bx];

    Lt_Zt_mats_px = Lt_Zt_mats[,proc_obs_ax:proc_obs_bx];

    if (!is.element("matrix", class(Lt_Zt_mats) ) || (nrow(Lt_Zt_mats) == 1) ) {
      Lt_Zt_mats_px = matrix(Lt_Zt_mats_px, 1);
    }

    y_tilde_px = y_tilde[proc_obs_ax:proc_obs_bx];
    neg_ll_px = 0;

    for (ix in proc_ax:proc_bx) {
      ax = grps_ixs_table[1,ix] - proc_obs_ax + 1;
      bx = grps_ixs_table[2,ix] - proc_obs_ax + 1;

      Lt_Zt_mats_ix = Lt_Zt_mats_px[,ax:bx];

      if (nrow(Lt_Zt_mats_px) == 1) {
        Lt_Zt_mats_ix = matrix(Lt_Zt_mats_ix, 1);
      }

      sigma_y_ix = tau2 * (crossprod(Lt_Zt_mats_ix) + diag(1, bx - ax + 1) );
      sigma_y_inv_ix = chol2inv(chol(sigma_y_ix) );

      neg_ll_px = (
        neg_ll_px +
        drop(
          crossprod(
            y_tilde_px[ax:bx],
            sigma_y_inv_ix %*% y_tilde_px[ax:bx]) ) +
        as.numeric(determinant(sigma_y_ix, logarithm=TRUE)$modulus) );
    }

    neg_ll_px
  };

  0.5 * ((m_dim * log(2 * pi) ) + neg_ll)
};


tucker_product_matrix = function(Zt_mats, L_list) {
  #   Description:
  #
  #     Arguments:
  #       Zt_mats: Matrix whose columns are the vectorized random effects covariate
  #                tensors.
  #
  #        L_list: List of D matrices, each of which is the square root covariance
  #                matrix to one of the D mode's random effects.
  #
  #  Return value:

  q_dims = unlist(lapply(L_list, nrow) );
  s_dims = unlist(lapply(L_list, ncol) );
  S_dim = prod(s_dims);

  Ks1q2_ixs = get_commutation_ixs(c(s_dims[1], q_dims[2]) );
  Ks2s1_ixs = get_commutation_ixs(rev(s_dims) );

  # This is equal to the Kronecker product of [L2t (x) L1t] times Zt_mats.
  Lt_Zt_mats = matrix(
    crossprod(
      L_list[[2]],
      matrix(
        matrix(
          crossprod(L_list[[1]], matrix(Zt_mats, q_dims[1]) ),
          q_dims[2] * s_dims[1])[Ks1q2_ixs,],
        q_dims[2]) ),
    S_dim);

  if (S_dim > 1) {
    Lt_Zt_mats = Lt_Zt_mats[Ks2s1_ixs,];
  }

  Lt_Zt_mats
};
