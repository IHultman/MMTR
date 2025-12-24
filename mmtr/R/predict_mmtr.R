#' This function computes predictions for an MMTR model.
#'
#' @param mmtr_mod MMTR model object produced by calling the mmtr() function.
#' @param X (D+1)-mode fixed effects tensor of dimension P1xP2x...xPDxM, i.e. there is a D-mode
#'   P1xP2x...xPD fixed effects covariate tensor associated with each of the M responses.
#' @param Z (D+1)-mode random effects tensor of dimension Q1xQ2x...xQDxM, i.e. there is a D-mode
#'   Q1xQ2x...xQD random effects covariate tensor associated with each of the M responses.
#' @param group_ids Integer or factor M-vector indicating which group each observation belongs to.
#' @export
predict_mmtr = function(
  mmtr_mod,
  X, Z,
  group_ids)
{
  #   Description:
  #
  #     Arguments:
  #        mmtr_mod: MMTR model object produced by calling the mmtr() function.
  #
  #               X: (D+1)-mode fixed effects tensor of dimension P1xP2x...xPDxM, i.e.
  #                  there is a D-mode P1xP2x...xPD fixed effects covariate tensor
  #                  associated with each of the M responses.
  #
  #               Z: (D+1)-mode random effects tensor of dimension Q1xQ2x...xQDxM, i.e.
  #                  there is a D-mode Q1xQ2x...xQD random effects covariate tensor
  #                  associated with each of the M responses.
  #
  #       group_ids: Integer or factor M-vector indicating which group each observation
  #                  belongs to.

  stopifnot(typeof(group_ids) == "integer");

  group_ids_test = as.factor(as.vector(group_ids) );
  all_group_levels = union(levels(group_ids_test), levels(mmtr_mod$group_ids) );
  group_ids_test = factor(group_ids_test, levels=all_group_levels);
  group_ids_train = factor(mmtr_mod$group_ids, levels=all_group_levels);

  stopifnot(is.element("array", class(X) ));
  stopifnot(is.element("array", class(Z) ));

  m_dim = length(group_ids_test);
  n_dim = length(unique(group_ids_test) );
  t_dim = length(dim(X) ) - 1;

  stopifnot(t_dim == 2);
  stopifnot(length(dim(Z) ) == (t_dim + 1) );
  stopifnot(dim(X)[t_dim+1] == m_dim);
  stopifnot(dim(Z)[t_dim+1] == m_dim);

  p_dims = dim(X)[1:t_dim];
  P_dim = prod(p_dims);
  q_dims = dim(Z)[1:t_dim];
  Q_dim = prod(q_dims);

  # Sort observations here so that group IDs are clustered. Save original ordering.
  grps_order_ixs = order(group_ids_test);
  grps_original_order_ixs = order(grps_order_ixs);
  group_ids_test = group_ids_test[grps_order_ixs];
  X = X[,,grps_order_ixs];
  Z = Z[,,grps_order_ixs];

  # Index ranges for each group.
  grps_info_table = as.data.frame(table(group_ids_test) );
  colnames(grps_info_table) = c("grp_id", "m_gx");
  grps_info_table$ax = cumsum(c(1, grps_info_table$m_gx[-n_dim]) );
  grps_info_table$bx = cumsum(grps_info_table$m_gx);

  # Check to make sure group ordering matches.
  for (gx in 1:n_dim) {
    grp_gx_ixs = grps_info_table$ax[gx]:grps_info_table$bx[gx];
    stopifnot(all(group_ids_test[grp_gx_ixs] == grps_info_table$grp_id[gx]) );
  }

  Xt_mats = matrix(X, P_dim, m_dim);
  Zt_mats = matrix(Z, Q_dim, m_dim);

  L_list = mmtr_mod$L_hat_list;

  stopifnot(typeof(L_list) == "list");
  stopifnot(length(L_list) == t_dim);

  for (kx in 1:t_dim) {
    stopifnot(is.element("matrix", class(L_list[[kx]]) ));
    stopifnot(nrow(L_list[[kx]]) == q_dims[kx]);
  }

  s_dims = unlist(lapply(L_list, ncol) );
  S_dim = prod(s_dims);

  Lt_Zt_mats = tucker_product_matrix(Zt_mats, L_list);

  y_hat = rep(NA, m_dim);

  for (gx in 1:n_dim) {
    gid_gx = grps_info_table$grp_id[gx];
    ax = grps_info_table$ax[gx];
    bx = grps_info_table$bx[gx];

    blup_gx = rep(0, S_dim);
    ix_blup = group_ids_train == gid_gx;

    if (any(ix_blup) ) {
      blup_gx = mmtr_mod$mu_hat_c_given_y[,ix_blup];
    }

    if (S_dim > 1) {
      y_hat[ax:bx] = crossprod(Lt_Zt_mats[,ax:bx], blup_gx);
    } else {
      y_hat[ax:bx] = Lt_Zt_mats[,ax:bx] * blup_gx;
    }
  }

  y_hat = y_hat + crossprod(Xt_mats, as.vector(mmtr_mod$B_hat) );

  y_hat[grps_original_order_ixs]
};
