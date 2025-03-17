library(stringr)
library(tools)

source("./JASA-QuasiLkld_lib.r");


args = commandArgs(TRUE);

if (length(args) < 1) {
  stop(
    paste(
      "No path to simulated data,",
      "simulation set index,",
      "or save directory provided.") );
}

if (length(args) < 3) {
  stop("Missing argument; three arguments expected.");
}

if (length(args) > 3) {
  stop(
    sprintf(
      paste0(
        "run_quasi_simulation_analysis.R takes exactly three arguments, ",
        "but %d arguments were provided."),
      length(args) ));
}

sim_data_dir = args[1];
ix_set = as.integer(args[2]);
save_dir = args[3];

if (!dir.exists(sim_data_dir) ) {
  stop(
    sprintf(
      "The provided data directory:\n%s\nwas not found.",
      sim_data_dir) );
}

if (!dir.exists(save_dir) ) {
  dir.create(save_dir, recursive=TRUE);

  stopifnot(dir.exists(save_dir) );
}

sim_sets = list.dirs(sim_data_dir, full.names=FALSE, recursive=FALSE);
n_sim_sets = length(sim_sets);

if ((ix_set < 1) || (ix_set > n_sim_sets) ) {
  stop(
    sprintf(
      "The provided file index %d is invalid. Valid arguments are integers between 1 and %d.",
      ix_set,
      n_sim_sets) );
}

sim_data_dir_sx = paste0(sim_data_dir, '/', sim_sets[ix_set], '/');
sim_filenames_sx = list.files(sim_data_dir_sx);
n_sim_filenames = length(sim_filenames_sx);

save_dir_sx = paste0(save_dir, '/', sim_sets[ix_set], '/');

if (!dir.exists(save_dir_sx) ) {
  dir.create(save_dir_sx);

  stopifnot(dir.exists(save_dir_sx) );
}

expected_fields = c(
  "X_arr",
  "Z_arr",
  "C_arr",
  "errs",
  "B_true",
  "L_true_list",
  "tau2_true",
  "y_vec",
  "group_ids");

a_scal_range = c(1e-4, 10);
n_a_vals = 10;
a_vals = exp(seq(log(a_scal_range[1]), log(a_scal_range[2]), length.out=n_a_vals) );

for (fx in 1:n_sim_filenames) {
  sim_data_full_path_fx = paste0(sim_data_dir_sx, sim_filenames_sx[fx]);

  if (file_ext(sim_data_full_path_fx) != "RData") {
    if (str_length(file_ext(sim_data_full_path_fx) ) > 0) {
      warning(
        sprintf(
          paste(
            "File of type RData required, ",
            "but filename with extension %s provided.",
            sep=""),
          file_ext(sim_data_full_path_fx) ));

    } else {
      warning("File of type RData required.");
    }

    next;
  }

  sim_data = readRDS(sim_data_full_path_fx);

  stopifnot(class(sim_data) == "list");

  stopifnot(
    (length(
      setdiff(
        expected_fields,
        names(sim_data) )) == 0) &&
    (length(
      setdiff(
        names(sim_data),
        expected_fields) ) == 0) );

  y = sim_data$y_vec;

  t_dim = 2;
  m_dim = length(y);

  stopifnot(length(dim(sim_data$X_arr) ) == (t_dim + 1) );
  stopifnot(length(dim(sim_data$Z_arr) ) == (t_dim + 1) );
  stopifnot(dim(sim_data$X_arr)[t_dim+1] == m_dim);
  stopifnot(dim(sim_data$Z_arr)[t_dim+1] == m_dim);

  P_dim = prod(dim(sim_data$X_arr)[1:t_dim]);
  Q_dim = prod(dim(sim_data$Z_arr)[1:t_dim]);

  X = t(matrix(sim_data$X_arr, P_dim, m_dim) );
  Z = t(matrix(sim_data$Z_arr, Q_dim, m_dim) );
  group_ids = sim_data$group_ids;

  set.seed(13);

  quasi_cv = CV_quasi_fixed_effects(X, y, Z, group_ids, a_vals);
  cv_test_mses = rowMeans(as.matrix(quasi_cv[,-1]) );
  a_best = max(quasi_cv$a[cv_test_mses == min(cv_test_mses)]);

  quasi_mod = Fixed_effects_estimation(X, y, Z, group_ids, a_best);

  save_filename = paste0(
    save_dir_sx,
    file_path_sans_ext(sim_filenames_sx[fx]),
    "_quasi_cv_results.RData");

  saveRDS(
    list(b_hat=as.numeric(quasi_mod$beta.hat) ),
    save_filename);
}
