library(dplyr)
library(ggplot2)
library(Matrix)
library(mmtr)
library(stringr)


sims_data_dir = paste0(
  "/home",
  "/ikhultman",
  "/Desktop",
  "/mmtr_paper_final",
  "/mmtr_paper_final",
  "/mmtr_analysis",
  "/mmtr_data",
  "/mmtr_equicorr_simulated_data/");

sims_results_dir = paste0(
  "/home",
  "/ikhultman",
  "/Desktop",
  "/mmtr_paper_final",
  "/mmtr_paper_final",
  "/mmtr_analysis",
  "/mmtr_results",
  "/mmtr_equicorr_simulated_data_results/");

sim_sets_dirs = list.files(sims_results_dir);
n_sim_sets = length(sim_sets_dirs);

sims_meta_data_filename = "simulations_meta_data.csv";
sims_meta_data = read.table(
  paste0(sims_data_dir, sims_meta_data_filename),
  header=TRUE,
  sep=",");

sim_results_tables = list();
alphas = rep(NA, n_sim_sets);

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
    data_sx_fx = readRDS(paste0(data_dir_sx, sim_results_meta_data$data_filename[fx]) );
    results_sx_fx = readRDS(paste0(results_dir_sx, sim_results_meta_data$results_filename[fx]) );

    if (fx == 1) {
      alphas[sx] = tcrossprod(data_sx_fx$L_gx)[2];
    }

    mmtr_mod_bic_min = results_sx_fx$mmtr_mod_min_ebic;

    sim_results_tables[[sx]]$B_f_norm[fx] = (
      norm(data_sx_fx$B_true - mmtr_mod_bic_min$B_hat, "F") /
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

    Zt_mats = matrix(data_sx_fx$X_arr, ncol=m_dim);
    L_hat = kronecker(mmtr_mod_bic_min$L_hat_list[[2]], mmtr_mod_bic_min$L_hat_list[[1]]);

    Lambda_hat_blocks = list();

    Lt_hat_Zt_mats = crossprod(L_hat, Zt_mats);

    for (gx in 1:n_dim) {
      ax = grps_ixs_table[1,gx];
      bx = grps_ixs_table[2,gx];
      m_gx = bx - ax + 1;

      Lt_hat_Zt_gx = matrix(Lt_hat_Zt_mats[,ax:bx], ncol=m_gx);
      Lambda_hat_blocks[[gx]] = mmtr_mod_bic_min$tau2_hat * (crossprod(Lt_hat_Zt_gx) + diag(1, m_gx) );
    }

    Lambda_hat = bdiag(Lambda_hat_blocks);
    Lambda_true = kronecker(diag(1, n_dim), tcrossprod(data_sx_fx$L_gx) );

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

alpha_order_ixs = order(alphas);
sim_norms_array = sim_norms_array[,,alpha_order_ixs];
alphas = alphas[alpha_order_ixs];

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

# Plots to see how the errors in estimating B change with alpha.
#plot(alphas, sim_norms_array[1,1,]);
#plot(alphas, sim_norms_array[2,1,]);
#plot(alphas, sim_norms_array[3,1,]);
#plot(alphas, sim_norms_array[4,1,]);
#plot(alphas, sim_norms_array[5,1,]);
#plot(alphas, sim_norms_array[6,1,]);
#plot(alphas, sim_norms_array[7,1,]);
#plot(alphas, sim_norms_array[8,1,]);

# Plots to see how the errors in estimating Sigma change with alpha.
#plot(alphas, sim_norms_array[1,2,]);
#plot(alphas, sim_norms_array[2,2,]);
#plot(alphas, sim_norms_array[3,2,]);
#plot(alphas, sim_norms_array[4,2,]);
#plot(alphas, sim_norms_array[5,2,]);
#plot(alphas, sim_norms_array[6,2,]);
#plot(alphas, sim_norms_array[7,2,]);
#plot(alphas, sim_norms_array[8,2,]);

#cx = 1;
#sims_meta_data_cx = sims_meta_data[sims_meta_data$caseID == cx,];
#n_sim_scenarios_cx = nrow(sims_meta_data_cx);
#plot_groups = sapply(
#  1:n_sim_scenarios_cx,
#  \(px) sprintf(
#    "Case %d, n = %d",
#    sims_meta_data_cx$caseID[px],
#    sims_meta_data_cx$n_dim[px]) );

#plot_df = data.frame(
#  sim_id=factor(rep(plot_groups, each=n_sim_sets), levels=plot_groups),
#  alpha=rep(alphas, n_sim_scenarios_cx),
#  lambda_err=as.vector(t(sim_norms_array[1:4,2,]) ));

#misspec_plots_diff_scales = (
#  ggplot(
#    data=plot_df,
#    mapping=aes(x=alpha, y=lambda_err) ) +
#  facet_wrap("sim_id", nrow=1, scales="free_y") +
#  geom_point() +
#  labs(
#    x="Alpha",
#    y="Lambda Relative Error") );

#misspec_plots_fixed_scales = (
#  ggplot(
#    data=plot_df,
#    mapping=aes(x=alpha, y=lambda_err) ) +
#  facet_wrap("sim_id", nrow=1) +
#  geom_point() +
#  labs(
#    x="Alpha",
#    y="Lambda Relative Error") );

plot_groups = sapply(
  1:n_sim_scenarios,
  \(px) sprintf(
    "Case %d, n = %d",
    sims_meta_data$caseID[px],
    sims_meta_data$n_dim[px]) );

plot_df = data.frame(
  sim_id=factor(rep(plot_groups, each=n_sim_sets), levels=plot_groups),
  alpha=rep(alphas, n_sim_scenarios),
  lambda_err=as.vector(t(sim_norms_array[,2,]) ));

misspec_plots_diff_scales = (
  ggplot(
    data=plot_df,
    mapping=aes(x=alpha, y=lambda_err) ) +
  facet_wrap("sim_id", nrow=2, scales="free_y") +
  geom_point() +
  labs(
    x="Alpha",
    y="Lambda Relative Error") );

diff_scales_save_filename = paste0(
  "/home",
  "/ikhultman",
  "/Desktop",
  "/mmtr_paper_final",
  "/mmtr_paper_final",
  "/results",
  "/misspec_diff_y_scales.png");

ggsave(
  diff_scales_save_filename,
  misspec_plots_diff_scales,
  width=18,
  height=12);

misspec_plots_fixed_scales = (
  ggplot(
    data=plot_df,
    mapping=aes(x=alpha, y=lambda_err) ) +
  facet_wrap("sim_id", nrow=2) +
  geom_point() +
  labs(
    x="Alpha",
    y="Lambda Relative Error") );

fixed_scales_save_filename = paste0(
  "/home",
  "/ikhultman",
  "/Desktop",
  "/mmtr_paper_final",
  "/mmtr_paper_final",
  "/results",
  "/misspec_fixed_y_scales.png");

ggsave(
  fixed_scales_save_filename,
  misspec_plots_fixed_scales,
  width=18,
  height=12);
