library(mmtr)
library(stringr)
library(tools)


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
        "run_mmtr_equicorr_simulation_analysis.R takes exactly three arguments, ",
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
  "B_true",
  "L_gx",
  "y_vec",
  "group_ids");

n_reps = 1000;
conv_thrsh = 0.001;

B_lambda_range = c(1e-4, 0.04);
L_lambda_range = c(1e-4, 1);

n_B_lambdas = 10;
n_L_lambdas = 10;

n_proc = n_B_lambdas;
s_dims_max = c(5, 5);

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
  X = sim_data$X_arr;
  Z = sim_data$X_arr;
  group_ids = sim_data$group_ids;

  set.seed(13);

  # Estimate model parameters over a grid of lambda values.
  lam_sweep_results = lambda_sweep_mmtr(
    y, X, Z,
    group_ids,
    s_dims_max,
    n_reps=n_reps,
    B_lambda_range=B_lambda_range,
    L_lambda_range=L_lambda_range,
    n_B_lambdas=n_B_lambdas,
    n_L_lambdas=n_L_lambdas,
    comp_ebic=TRUE,
    conv_thrsh=conv_thrsh,
    n_proc=n_proc);

  save_filename = paste0(
    save_dir_sx,
    file_path_sans_ext(sim_filenames_sx[fx]),
    "_lambda_sweep_results.RData");

  saveRDS(lam_sweep_results, save_filename);
}
