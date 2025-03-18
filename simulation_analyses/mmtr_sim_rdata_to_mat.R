library(R.matlab)
library(tools)


data_dir = "./simulated_mmtr_data/";
save_dir = "./simulated_mmtr_matlab_data/";

if (!dir.exists(save_dir) ) {
  dir.create(save_dir, recursive=TRUE);

  stopifnot(dir.exists(save_dir) );
}

sim_dirs = list.dirs(data_dir, full.name=FALSE, recursive=FALSE);
n_sim_dirs = length(sim_dirs);

for (dx in 1:n_sim_dirs) {
  next_sim_dir = sim_dirs[dx];
  next_save_dir = paste0(save_dir, "/", next_sim_dir);

  if (!dir.exists(next_save_dir) ) {
    dir.create(next_save_dir);

    stopifnot(dir.exists(next_save_dir) );
  }

  next_sim_fileset = list.files(paste0(data_dir, "/", next_sim_dir) );
  n_sim_files = length(next_sim_fileset);

  for (fx in 1:n_sim_files) {
    next_sim_file = next_sim_fileset[fx];
    next_sim_data = readRDS(
      paste(data_dir, next_sim_dir, next_sim_file, sep="/") );

    next_save_filename = paste0(
      next_save_dir, "/", file_path_sans_ext(next_sim_file), ".mat");

    writeMat(
      next_save_filename,
      X=next_sim_data$X_arr,
      Z=next_sim_data$Z_arr,
      y=next_sim_data$y_vec,
      group_ids=next_sim_data$group_ids,
      B=next_sim_data$B_true,
      Sigma_1=tcrossprod(next_sim_data$L_true_list[[1]]),
      Sigma_2=tcrossprod(next_sim_data$L_true_list[[2]]),
      tau2=next_sim_data$tau2_true);
  }
}
