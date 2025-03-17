set.seed(13);

# Set the number of modes for fixed and random effects covariate tensors.
t_dim = 2;

# Set the different parameter tensor dimensions corresponding to each
# simulation case.
p_dims = matrix(
  c(c(5, 5),
    c(10, 10) ),
  ncol=t_dim,
  byrow=TRUE);

n_cases = nrow(p_dims);

# Set the different number of groups for which to simulate data for each case.
n_dim = t(matrix(
  c(18, 27, 54, 81, 64, 96, 192, 288),
  ncol=n_cases) );

# Set the number of observations to simulate for each different group size
# for each case.
n_grp_obs = 6;

n_sim_sizes = ncol(n_dim);
n_obs_sizes = length(n_grp_obs);

cases_df = data.frame(cbind(1:n_cases, p_dims) );
colnames(cases_df) = c(
  "caseID",
  paste("p_dim", 1:t_dim, sep="") );

groups_df = data.frame(
  group_sizeID=rep(1:n_sim_sizes, each=n_cases),
  caseID=rep(1:n_cases, n_sim_sizes),
  n_dim=as.vector(n_dim) );

obs_df = data.frame(
  obs_sizeID=rep(1:n_obs_sizes, each=n_cases),
  caseID=rep(1:n_cases, n_obs_sizes),
  n_grp_obs=n_grp_obs);

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

save_dir = paste0("./simulated_equicorr_data/");

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

  # Generate random equicorrelation alpha value.
  equicorr_alpha = runif(1, 0.2, 0.8);

  L_true_list = list();

  for (osx in 1:n_obs_sizes) {
    n_obs_per_group = n_grp_obs[osx];
    sigma_gx = matrix(equicorr_alpha, n_obs_per_group, n_obs_per_group);
    diag(sigma_gx) = 1;

    L_true_list[[osx]] = t(chol(sigma_gx) );
  }

  X_list = list();

  group_sizes_list =  list();
  group_ixs_list = list();

  for (cx in 1:n_cases) {
    X_list[[cx]] = list();

    group_sizes_list[[cx]] = list();
    group_ixs_list[[cx]] = list();

    p_dims_local = p_dims[cx,];
    P_dim_local = prod(p_dims_local);

    for (gsx in 1:n_sim_sizes) {
      X_list[[cx]][[gsx]] = list();

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

        obs_per_grp_local = sim_pars_table$n_grp_obs[m_dim_ix];
        group_sizes_list[[cx]][[gsx]][[osx]] = rep(obs_per_grp_local, n_dim_local);

        group_ixs_list[[cx]][[gsx]][[osx]] = cbind(
          cumsum(c(1, group_sizes_list[[cx]][[gsx]][[osx]][-n_dim_local]) ),
          cumsum(group_sizes_list[[cx]][[gsx]][[osx]]) );

        stopifnot(m_dim_local == (obs_per_grp_local * n_dim_local) );
        stopifnot(nrow(group_ixs_list[[cx]][[gsx]][[osx]]) == n_dim_local);
        stopifnot(sum(group_sizes_list[[cx]][[gsx]][[osx]]) == m_dim_local);
      }
    }
  }

  # Consistency and correctness checks for data simulated in the previous code block.
  for (cx in 1:n_cases) {
    for (gsx in 1:n_sim_sizes) {
      X_arr_sizes = unlist(
        lapply(
          X_list[[cx]][[gsx]],
          function(X_arr) dim(X_arr)[t_dim+1]) );

      grp_sz_totals_check = unlist(lapply(group_sizes_list[[cx]][[gsx]], sum) );

      m_dim_ixs = which(
        (sim_pars_table$caseID == cx) &
        (sim_pars_table$group_sizeID == gsx) );

      m_dims_check = sim_pars_table$m_dim[m_dim_ixs];

      stopifnot(all(X_arr_sizes == m_dims_check) );
      stopifnot(all(grp_sz_totals_check == m_dims_check) );

      grp_sz_vec_sizes = unlist(lapply(group_sizes_list[[cx]][[gsx]], length) );
      grp_ixs_check = unlist(lapply(group_ixs_list[[cx]][[gsx]], nrow) );

      n_dim_ix = which(
        (groups_df$caseID == cx) &
        (groups_df$group_sizeID == gsx) );

      n_dim_check = groups_df$n_dim[n_dim_ix];

      stopifnot(all(grp_sz_vec_sizes == n_dim_check) );
      stopifnot(all(grp_ixs_check == n_dim_check) );
    }
  }

  # Simulate responses.
  y_list = list();

  for (cx in 1:n_cases) {
    y_list[[cx]] = list();

    p_dims_local = p_dims[cx,];
    P_dim_local = prod(p_dims_local);

    for (gsx in 1:n_sim_sizes) {
      y_list[[cx]][[gsx]] = list();

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

        errs = L_true_list[[osx]] %*% matrix(rnorm(m_dim_local), ncol=n_dim_local);
        X_mat = t(matrix(X_list[[cx]][[gsx]][[osx]], P_dim_local, m_dim_local) );
        y_list[[cx]][[gsx]][[osx]] = (X_mat %*% as.vector(B_true_list[[cx]]) ) + as.vector(errs);
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
          B_true=B_true_list[[cx]],
          L_gx=L_true_list[[osx]],
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
