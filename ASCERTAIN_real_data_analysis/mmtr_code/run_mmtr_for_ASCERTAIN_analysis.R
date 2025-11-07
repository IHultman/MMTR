library(mmtr)
library(stringr)
library(tools)


args = commandArgs(TRUE);

if (length(args) < 1) {
  stop(
    paste(
      "No path to ASCERTAIN data,",
      "file index,",
      "response variable label",
      "or save directory provided.") );
}

if (length(args) < 4) {
  stop("Missing argument; four arguments expected.");
}

if (length(args) > 4) {
  stop(
    sprintf(
      paste0(
        "run_mmtr_for_ASCERTAIN_analysis.R takes exactly four arguments, ",
        "but %d arguments were provided."),
      length(args) ));
}

ASCERTAIN_data_dir = args[1];
file_ix = as.integer(args[2]);
ASCERTAIN_score_ix = as.integer(args[3]);
save_dir = args[4];

if (!dir.exists(ASCERTAIN_data_dir) ) {
  stop(
    sprintf(
      "The provided data directory %s was not found.",
      ASCERTAIN_data_dir) );
}

ASCERTAIN_data_filenames = list.files(ASCERTAIN_data_dir);
n_ASCERTAIN_data_files = length(ASCERTAIN_data_filenames);

if ((file_ix < 1) || (file_ix > n_ASCERTAIN_data_files) ) {
  stop(
    sprintf(
      "The provided file index %d is invalid. Valid arguments are integers between 1 and %d.",
      file_ix,
      n_ASCERTAIN_data_files) );
}

ASCERTAIN_data_filename_fx = ASCERTAIN_data_filenames[file_ix];
ASCERTAIN_data_full_path_fx = file.path(ASCERTAIN_data_dir, ASCERTAIN_data_filename_fx);

if (file_ext(ASCERTAIN_data_filename_fx) != "RData") {
  if (str_length(file_ext(ASCERTAIN_data_filename_fx) ) > 0) {
    stop(
      sprintf(
        paste(
          "File of type RData required, ",
          "but filename with extension %s provided.",
          sep=""),
        file_ext(ASCERTAIN_data_filename_fx) ));

  } else {
    stop("File of type RData required.");
  }
}

if (!file.exists(save_dir) ) {
  dir.create(save_dir, recursive=TRUE);
}

stopifnot(dir.exists(save_dir) );

ASCERTAIN_data = readRDS(ASCERTAIN_data_full_path_fx);

stopifnot(class(ASCERTAIN_data) == "list");

expected_fields = c(
  "X_train",
  "X_test",
  "y_train",
  "y_test",
  "group_ids_train",
  "group_ids_test");

stopifnot(
  (length(
    setdiff(
      expected_fields,
      names(ASCERTAIN_data) )) == 0) &&
  (length(
    setdiff(
      names(ASCERTAIN_data),
      expected_fields) ) == 0) );

if ((ASCERTAIN_score_ix < 1) || (ASCERTAIN_score_ix > ncol(ASCERTAIN_data$y_train) )) {
  stop(
    sprintf(
      "The provided response matrix column index %d is invalid. Valid arguments are integers between 1 and %d.",
      ASCERTAIN_score_ix,
      ncol(ASCERTAIN_data$y_train) ));
}

y = ASCERTAIN_data$y_train[,ASCERTAIN_score_ix];
X = ASCERTAIN_data$X_train;
Z = ASCERTAIN_data$X_train;
group_ids = ASCERTAIN_data$group_ids_train;

n_reps = 1000;
conv_thrsh = 0.001;

B_lambda_range = c(1e-4, 0.05);
L_lambda_range = c(1e-4, 1);

n_B_lambdas = 10;
n_L_lambdas = 10;

n_proc = 10;
s_dims_max = c(3, 3);

set.seed(13);

# Estimate model parameters over a grid of lambda values.
ASCERTAIN_lam_sweep_results = lambda_sweep_mmtr(
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

save_filename = file.path(
  save_dir,
  sprintf(
    "%s_y_ix_%d_lambda_sweep_results.RData",
    file_path_sans_ext(ASCERTAIN_data_filename_fx),
    ASCERTAIN_score_ix) );

saveRDS(ASCERTAIN_lam_sweep_results, save_filename);
