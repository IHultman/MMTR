set.seed(13);

# Set the number of modes for fixed and random effects covariate tensors.
t_dim = 2;

# Set the different parameter tensor dimensions corresponding to each
# simulation case. The `p_dims`, `q_dims` and `s_dims` matrices should
# have the same number of rows and columns. The number of columns must be
# equal to the `t_dim` variable set above.
p_dims = matrix(
  c(c(5, 5),
    c(10, 10) ),
  ncol=t_dim,
  byrow=TRUE);

q_dims =  p_dims;

s_dims = matrix(
  c(c(2, 2),
    c(3, 3) ),
  ncol=t_dim,
  byrow=TRUE);

n_cases = nrow(p_dims);

stopifnot(nrow(q_dims) == n_cases);
stopifnot(nrow(s_dims) == n_cases);

# Set the different number of groups for which to simulate data for each case.
n_case_params_vec = apply(p_dims, 1, prod) + rowSums(q_dims * s_dims);
percent_n_dim = c(0.4, 0.6, 1.2, 1.8);
n_dim = round(
  matrix(rep(percent_n_dim, each=n_cases), n_cases) *
  matrix(rep(n_case_params_vec, length(percent_n_dim) ), n_cases) );

# Set the number of observations to simulate for each different group size
# for each case.
n_grp_obs = matrix(
  rep(c(2, 6), n_cases),
  nrow=n_cases,
  byrow=TRUE); # Number of observations in each group.

n_sim_sizes = ncol(n_dim);
n_obs_sizes = ncol(n_grp_obs);

cases_df = data.frame(cbind(1:n_cases, p_dims, q_dims, s_dims) );
colnames(cases_df) = c(
  "caseID",
  paste("p_dim", 1:t_dim, sep=""),
  paste("q_dim", 1:t_dim, sep=""),
  paste("s_dim", 1:t_dim, sep="") );

groups_df = data.frame(
  group_sizeID=rep(1:n_sim_sizes, each=n_cases),
  caseID=rep(1:n_cases, n_sim_sizes),
  n_dim=as.vector(n_dim) );

obs_df = data.frame(
  obs_sizeID=rep(1:n_obs_sizes, each=n_cases),
  caseID=rep(1:n_cases, n_obs_sizes),
  n_grp_obs=as.vector(n_grp_obs) );

# Each row of this table gives a different combination of parameters to use
# for different simulations.
sim_pars_table = expand.grid(
  group_sizeID=unique(groups_df$group_sizeID),
  obs_sizeID=unique(obs_df$obs_sizeID) );

sim_pars_table = merge(
  merge(sim_pars_table, obs_df),
  groups_df);

sim_pars_table$m_dim = (
  sim_pars_table$n_grp_obs * sim_pars_table$n_dim);

sim_pars_table = merge(sim_pars_table, cases_df);

sim_pars_ixs = order(
  sim_pars_table$caseID,
  sim_pars_table$group_sizeID,
  sim_pars_table$obs_sizeID);

sim_pars_table = sim_pars_table[sim_pars_ixs,];
rownames(sim_pars_table) = NULL;

save_dir = paste0("./simulated_mmtr_data/");

if (!dir.exists(save_dir) ) {
  dir.create(save_dir, recursive=TRUE);

  stopifnot(dir.exists(save_dir) );
}

write.table(
  sim_pars_table,
  file=paste0(save_dir, "/simulations_meta_data.csv"),
  row.names=FALSE,
  sep=", ");

n_sim_sets = 100;

for (sx in 1:n_sim_sets) {
  # Set true B parameter array values for each case.
  B_true_list = list();
  percent_nonzero = 0.60;

  for (cx in 1:n_cases) {
    B_true_list[[cx]] = array(0, dim=p_dims[cx,]);

    n_arr_elems = prod(p_dims[cx,]);
    n_nonzero = floor(percent_nonzero * n_arr_elems);
    B_nz_ixs = sort(sample(1:n_arr_elems, n_nonzero) );

    B_true_list[[cx]][B_nz_ixs] = sample(
      seq(-10, 10, 0.5),
      n_nonzero,
      replace=TRUE);
  }

  # Set true sigma covariance matrix values for each mode and each case
  # (saved as square root matrices).
  L_true_list = list();

  for (cx in 1:n_cases) {
    L_true_list[[cx]] = list();

    for (dx in 1:t_dim) {
      sig_rnk = s_dims[cx,dx];

      L_cx_dx = matrix(
        runif(q_dims[cx,dx] * sig_rnk, -1, 1),
        q_dims[cx,dx],
        sig_rnk);

      #L_cx_dx_svd = svd(L_cx_dx);
      #diag_scale_val = exp((1 / sig_rnk) * sum(log(L_cx_dx_svd$d) ));
      #L_true_list[[cx]][[dx]] = L_cx_dx_svd$u %*% diag(L_cx_dx_svd$d, sig_rnk) / diag_scale_val;

      # Create random variance/covariance matrix for the dth mode with rank s_d.
      mode_d_sigma = tcrossprod(L_cx_dx);

      # Scale the variance/covariance matrix such that the product of the square
      # roots of its nonzero singular values is one.
      mode_d_sig_svd = svd(mode_d_sigma);
      diag_scale_val = exp((1 / sig_rnk) * sum(log(mode_d_sig_svd$d[1:sig_rnk]) ));
      sig_sqrt_diag = diag(sqrt(mode_d_sig_svd$d[1:sig_rnk] / diag_scale_val) );
      L_true_list[[cx]][[dx]] = mode_d_sig_svd$u[,1:sig_rnk] %*% sig_sqrt_diag;
    }
  }

  # Set true tau2 parameter value.
  tau2_true = 0.5;
  tau_true = sqrt(tau2_true);

  # For each case and each group size, construct the fixed and random
  # effects covariate arrays along with the observed random effects
  # arrays.
  X_list = list();
  C_list = list();

  errs_list = list();
  group_sizes_list =  list();
  group_ixs_list = list();

  for (cx in 1:n_cases) {
    X_list[[cx]] = list();
    C_list[[cx]] = list();

    errs_list[[cx]] = list();
    group_sizes_list[[cx]] = list();
    group_ixs_list[[cx]] = list();

    p_dims_local = p_dims[cx,];
    q_dims_local = q_dims[cx,];
    s_dims_local = s_dims[cx,];

    P_dim_local = prod(p_dims_local);
    Q_dim_local = prod(q_dims_local);
    S_dim_local = prod(s_dims_local);

    for (gsx in 1:n_sim_sizes) {
      X_list[[cx]][[gsx]] = list();

      errs_list[[cx]][[gsx]] = list();
      group_sizes_list[[cx]][[gsx]] = list();
      group_ixs_list[[cx]][[gsx]] = list();

      n_dim_ix = which(
        (groups_df$caseID == cx) &
        (groups_df$group_sizeID == gsx) );

      n_dim_local = groups_df$n_dim[n_dim_ix];

      for (osx in 1:n_obs_sizes) {
        m_dim_ix = which(
          (sim_pars_table$caseID == cx) &
          (sim_pars_table$group_sizeID == gsx) &
          (sim_pars_table$obs_sizeID == osx) );

        m_dim_local = sim_pars_table$m_dim[m_dim_ix];

        X_list[[cx]][[gsx]][[osx]] = array(
          rnorm(P_dim_local * m_dim_local),
          dim=c(p_dims_local, m_dim_local) );

        errs_list[[cx]][[gsx]][[osx]] = rnorm(m_dim_local, sd=tau_true);

        obs_per_grp_local = sim_pars_table$n_grp_obs[m_dim_ix];
        group_sizes_list[[cx]][[gsx]][[osx]] = rep(obs_per_grp_local, n_dim_local);

        group_ixs_list[[cx]][[gsx]][[osx]] = cbind(
          cumsum(c(1, group_sizes_list[[cx]][[gsx]][[osx]][-n_dim_local]) ),
          cumsum(group_sizes_list[[cx]][[gsx]][[osx]]) );

        stopifnot(m_dim_local == (obs_per_grp_local * n_dim_local) );
        stopifnot(nrow(group_ixs_list[[cx]][[gsx]][[osx]]) == n_dim_local);
        stopifnot(sum(group_sizes_list[[cx]][[gsx]][[osx]]) == m_dim_local);
      }

      C_list[[cx]][[gsx]] = matrix(
        rnorm(S_dim_local * n_dim_local, sd=tau_true),
        S_dim_local,
        n_dim_local);
    }
  }

  Z_list = X_list;

  # Consistency and correctness checks for data simulated in the previous code block.
  for (cx in 1:n_cases) {
    for (gsx in 1:n_sim_sizes) {
      X_arr_sizes = unlist(
        lapply(
          X_list[[cx]][[gsx]],
          function(X_arr) dim(X_arr)[t_dim+1]) );

      Z_arr_sizes = unlist(
        lapply(
          Z_list[[cx]][[gsx]],
          function(Z_arr) dim(Z_arr)[t_dim+1]) );

      err_vec_sizes = unlist(lapply(errs_list[[cx]][[gsx]], length) );
      grp_sz_totals_check = unlist(lapply(group_sizes_list[[cx]][[gsx]], sum) );

      m_dim_ixs = which(
        (sim_pars_table$caseID == cx) &
        (sim_pars_table$group_sizeID == gsx) );

      m_dims_check = sim_pars_table$m_dim[m_dim_ixs];

      stopifnot(all(X_arr_sizes == m_dims_check) );
      stopifnot(all(Z_arr_sizes == m_dims_check) );
      stopifnot(all(err_vec_sizes == m_dims_check) );
      stopifnot(all(grp_sz_totals_check == m_dims_check) );

      grp_sz_vec_sizes = unlist(lapply(group_sizes_list[[cx]][[gsx]], length) );
      grp_ixs_check = unlist(lapply(group_ixs_list[[cx]][[gsx]], nrow) );

      n_dim_ix = which(
        (groups_df$caseID == cx) &
        (groups_df$group_sizeID == gsx) );

      n_dim_check = groups_df$n_dim[n_dim_ix];

      stopifnot(all(grp_sz_vec_sizes == n_dim_check) );
      stopifnot(all(grp_ixs_check == n_dim_check) );
      stopifnot(ncol(C_list[[cx]][[gsx]]) == n_dim_check);
    }
  }

  # Simulate responses.
  y_list = list();

  for (cx in 1:n_cases) {
    y_list[[cx]] = list();

    p_dims_local = p_dims[cx,];
    q_dims_local = q_dims[cx,];

    P_dim_local = prod(p_dims_local);
    Q_dim_local = prod(q_dims_local);

    L_mat = 1;

    for (kx in t_dim:1) {
      L_mat = kronecker(L_mat, L_true_list[[cx]][[kx]]);
    }

    for (gsx in 1:n_sim_sizes) {
      y_list[[cx]][[gsx]] = list();

      n_dim_ix = which(
        (groups_df$caseID == cx) &
        (groups_df$group_sizeID == gsx) );

      n_dim_local = groups_df$n_dim[n_dim_ix];

      A_mat = L_mat %*% C_list[[cx]][[gsx]];

      for (osx in 1:n_obs_sizes) {
        m_dim_ix = which(
          (sim_pars_table$caseID == cx) &
          (sim_pars_table$group_sizeID == gsx) &
          (sim_pars_table$obs_sizeID == osx) );

        m_dim_local = sim_pars_table$m_dim[m_dim_ix];

        y_list[[cx]][[gsx]][[osx]] = rep(NA, m_dim_local);

        group_ixs_local = group_ixs_list[[cx]][[gsx]][[osx]];

        for (grp_ix in 1:n_dim_local) {
          for (mx in group_ixs_local[grp_ix,1]:group_ixs_local[grp_ix,2]) {
            x_ax = ((mx - 1) * P_dim_local) + 1;
            x_bx = mx * P_dim_local;

            z_ax = ((mx - 1) * Q_dim_local) + 1;
            z_bx = mx * Q_dim_local;

            y_list[[cx]][[gsx]][[osx]][mx] = (
              sum(
                X_list[[cx]][[gsx]][[osx]][x_ax:x_bx] *
                as.vector(B_true_list[[cx]]) ) +
              sum(Z_list[[cx]][[gsx]][[osx]][z_ax:z_bx] * A_mat[,grp_ix]) +
              errs_list[[cx]][[gsx]][[osx]][mx]);
          }
        }
      }
    }
  }

  # Final consistency check.
  for (cx in 1:n_cases) {
    for (gsx in 1:n_sim_sizes) {
      y_sizes = unlist(lapply(y_list[[cx]][[gsx]], length) );

      m_dim_ixs = which(
        (sim_pars_table$caseID == cx) &
        (sim_pars_table$group_sizeID == gsx) );

      m_dims_check = sim_pars_table$m_dim[m_dim_ixs];

      stopifnot(all(y_sizes == m_dims_check) );
    }
  }

  set_save_dir = paste0(save_dir, '/sim_set_', sx);

  if (!dir.exists(set_save_dir) ) {
    dir.create(set_save_dir);

    stopifnot(dir.exists(set_save_dir) );
  }

  for (cx in 1:n_cases) {
    for (gsx in 1:n_sim_sizes) {
      for (osx in 1:n_obs_sizes) {
        next_group_ids = rep(
          1:length(group_sizes_list[[cx]][[gsx]][[osx]]),
          group_sizes_list[[cx]][[gsx]][[osx]]);

        sim_instance_pars = list(
          X_arr=X_list[[cx]][[gsx]][[osx]],
          Z_arr=Z_list[[cx]][[gsx]][[osx]],
          C_arr=C_list[[cx]][[gsx]],
          errs=errs_list[[cx]][[gsx]][[osx]],
          B_true=B_true_list[[cx]],
          L_true_list=L_true_list[[cx]],
          tau2_true=tau2_true,
          y_vec=y_list[[cx]][[gsx]][[osx]],
          group_ids=next_group_ids);

        save_filename = paste0(
          set_save_dir,
          sprintf(
            "/simulation_caseID_%d_group_sizeID_%d_obs_sizeID_%d.RData",
            cx, gsx, osx) );

        saveRDS(sim_instance_pars, save_filename);
      }
    }
  }
}
