library(dplyr)
library(Matrix)
library(R.matlab)
library(stringr)


sims_data_dir = paste0(
  "/home",
  "/ikhultman",
  "/Desktop",
  "/mmtr_paper_final",
  "/mmtr_paper_final",
  "/mmtr_analysis",
  "/mmtr_data",
  "/mmtr_simulated_data/");

sims_results_dir = paste0(
  "/home",
  "/ikhultman",
  "/Desktop",
  "/mmtr_paper_final",
  "/mmtr_paper_final",
  "/kruskal_analysis",
  "/kruskal_results",
  "/gee_simulated_data_results_equicorr/");

sim_sets_dirs = list.files(sims_results_dir);
n_sim_sets = length(sim_sets_dirs);

sims_meta_data_filename = "simulations_meta_data.csv";
sims_meta_data = read.table(
  paste0(sims_data_dir, sims_meta_data_filename),
  header=TRUE,
  sep=",");

sim_results_tables = list();

for (sx in 1:n_sim_sets) {
  data_dir_sx = paste0(sims_data_dir, '/', sim_sets_dirs[sx], '/');
  results_dir_sx = paste0(sims_results_dir, '/', sim_sets_dirs[sx], '/');

  sim_results_meta_data = merge(
    data.frame(results_filename=list.files(results_dir_sx) ) %>%
    mutate(
      caseID=as.numeric(str_extract(results_filename, "(?<=caseID_)[0-9]*") ),
      group_sizeID=as.numeric(str_extract(results_filename, "(?<=group_sizeID_)[0-9]*") ),
      obs_sizeID=as.numeric(str_extract(results_filename, "(?<=obs_sizeID_)[0-9]*") )),
    data.frame(data_filename=list.files(data_dir_sx) ) %>%
    mutate(
      caseID=as.numeric(str_extract(data_filename, "(?<=caseID_)[0-9]*") ),
      group_sizeID=as.numeric(str_extract(data_filename, "(?<=group_sizeID_)[0-9]*") ),
      obs_sizeID=as.numeric(str_extract(data_filename, "(?<=obs_sizeID_)[0-9]*") )));

  sim_results_meta_data = merge(sims_meta_data, sim_results_meta_data);

  sims_ord_ixs = order(
    sim_results_meta_data$caseID,
    sim_results_meta_data$group_sizeID,
    sim_results_meta_data$obs_sizeID);

  sim_results_meta_data = sim_results_meta_data[sims_ord_ixs,];
  rownames(sim_results_meta_data) = NULL;

  n_sims = nrow(sim_results_meta_data);

  sim_results_tables[[sx]] = data.frame(
    caseID=sim_results_meta_data$caseID,
    group_sizeID=sim_results_meta_data$group_sizeID,
    obs_sizeID=sim_results_meta_data$obs_sizeID,
    B_f_norm=double(n_sims),
    Lambda_f_norm=double(n_sims) );

  for (fx in 1:n_sims) {
    data_sx_fx = readRDS(
      paste0(data_dir_sx, sim_results_meta_data$data_filename[fx]) );

    results_sx_fx = readMat(
      paste0(results_dir_sx, sim_results_meta_data$results_filename[fx]) );

    p_dims = dim(data_sx_fx$B_true);
    B_hat = array(results_sx_fx$betahat, dim=p_dims);

    sim_results_tables[[sx]]$B_f_norm[fx] = (
      norm(data_sx_fx$B_true - B_hat, "F") /
      norm(data_sx_fx$B_true, "F") );

    m_dim = length(data_sx_fx$y_vec);
    n_dim = length(unique(data_sx_fx$group_ids) );

    # Index ranges for each group.
    grps_info_table = as.data.frame(table(data_sx_fx$group_ids) );
    colnames(grps_info_table) = c("grp_id", "m_gx");
    grps_info_table$ax = cumsum(c(1, grps_info_table$m_gx[-n_dim]) );
    grps_info_table$bx = cumsum(grps_info_table$m_gx);

    # Check to make sure group ordering matches.
    for (gx in 1:n_dim) {
      grp_gx_ixs = grps_info_table$ax[gx]:grps_info_table$bx[gx];
      stopifnot(all(data_sx_fx$group_ids[grp_gx_ixs] == grps_info_table$grp_id[gx]) );
    }

    grps_ixs_table = t(as.matrix(grps_info_table[,c("ax", "bx")]) );
    rownames(grps_ixs_table) = NULL;

    Zt_mats = matrix(data_sx_fx$Z_arr, ncol=m_dim);
    L_mat = kronecker(data_sx_fx$L_true_list[[2]], data_sx_fx$L_true_list[[1]]);
    Lambda_blocks = list();

    for (gx in 1:n_dim) {
      ax = grps_ixs_table[1,gx];
      bx = grps_ixs_table[2,gx];
      m_gx = bx - ax + 1;

      Lt_Zt_gx = crossprod(L_mat, Zt_mats[,ax:bx]);
      Lambda_blocks[[gx]] = data_sx_fx$tau2_true * (crossprod(Lt_Zt_gx) + diag(1, m_gx) );
    }

    Lambda_true = bdiag(Lambda_blocks);
    Lambda_hat = kronecker(Diagonal(n_dim), results_sx_fx$alphahat);

    sim_results_tables[[sx]]$Lambda_f_norm[fx] = (
      norm(Lambda_true - Lambda_hat, "F") /
      norm(Lambda_true, "F") );
  }
}

# Number of significant digits to round to.
n_signif = 2;

group_info_colnames = c("caseID", "group_sizeID", "obs_sizeID");
norm_colnames = c(
  "B_f_norm",
  "Lambda_f_norm");

order_check_mat = sim_results_tables[[1]][,group_info_colnames];

for (sx in 2:n_sim_sets) {
  stopifnot(
    identical(
      sim_results_tables[[sx]][,group_info_colnames],
      order_check_mat) );
}

n_norms = length(norm_colnames);
n_sim_scenarios = nrow(sims_meta_data);
sim_norms_array = array(NA, dim=c(n_sim_scenarios, n_norms, n_sim_sets) );

for (sx in 1:n_sim_sets) {
  sim_norms_array[,,sx] = as.matrix(sim_results_tables[[sx]][,norm_colnames]);
}

mean_sim_results_table = as.data.frame(
  signif(
    apply(sim_norms_array, MARGIN=c(1, 2), mean, na.rm=TRUE),
    n_signif) );

colnames(mean_sim_results_table) = norm_colnames;

mean_sim_results_table = cbind(
  sim_results_tables[[1]][,group_info_colnames],
  mean_sim_results_table);

sd_sim_results_table = as.data.frame(
  signif(
    apply(sim_norms_array, MARGIN=c(1, 2), sd, na.rm=TRUE),
    n_signif) );

colnames(sd_sim_results_table) = norm_colnames;

sd_sim_results_table = cbind(
  sim_results_tables[[1]][,group_info_colnames],
  sd_sim_results_table);
