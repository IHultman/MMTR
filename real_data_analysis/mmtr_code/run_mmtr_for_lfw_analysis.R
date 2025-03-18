library(mmtr)
library(stringr)
library(tools)


args = commandArgs(TRUE);

if (length(args) < 1) {
  stop(
    paste(
      "No path to lfw data,",
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
        "run_simulation_paper_final_lfw.R takes exactly four arguments, ",
        "but %d arguments were provided."),
      length(args) ));
}

lfw_data_dir = args[1];
file_ix = as.integer(args[2]);
lfw_attribute = args[3];
save_dir = args[4];

if (!dir.exists(lfw_data_dir) ) {
  stop(
    sprintf(
      "The provided data directory %s was not found.",
      lfw_data_dir) );
}

lfw_data_filenames = list.files(lfw_data_dir);
n_lfw_data_files = length(lfw_data_filenames);

if ((file_ix < 1) || (file_ix > n_lfw_data_files) ) {
  stop(
    sprintf(
      "The provided file index %d is invalid. Valid arguments are integers between 1 and %d.",
      file_ix,
      n_lfw_data_files) );
}

lfw_data_filename_fx = lfw_data_filenames[file_ix];
lfw_data_full_path_fx = paste0(lfw_data_dir, "/", lfw_data_filename_fx);

if (file_ext(lfw_data_filename_fx) != "RData") {
  if (str_length(file_ext(lfw_data_filename_fx) ) > 0) {
    stop(
      sprintf(
        paste(
          "File of type RData required, ",
          "but filename with extension %s provided.",
          sep=""),
        file_ext(lfw_data_filename_fx) ));

  } else {
    stop("File of type RData required.");
  }
}

lfw_data = readRDS(lfw_data_full_path_fx);

stopifnot(class(lfw_data) == "list");

expected_fields = c(
  "X_train",
  "X_test",
  "y_train",
  "y_test",
  "group_labels_train",
  "group_labels_test");

stopifnot(
  (length(
    setdiff(
      expected_fields,
      names(lfw_data) )) == 0) &&
  (length(
    setdiff(
      names(lfw_data),
      expected_fields) ) == 0) );

stopifnot(is.element(lfw_attribute, colnames(lfw_data$y_train) ));

y = lfw_data$y_train[[lfw_attribute]];
X = lfw_data$X_train;
Z = lfw_data$X_train;
group_ids = lfw_data$group_labels_train;

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
lfw_lam_sweep_results = lambda_sweep_mmtr(
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

if (!dir.exists(save_dir) ) {
  dir.create(save_dir, recursive=TRUE);

  stopifnot(dir.exists(save_dir) );
}

save_filename = paste0(
  save_dir,
  sprintf(
    "/%s_%s_lambda_sweep_results.RData",
    file_path_sans_ext(lfw_data_filename_fx),
    lfw_attribute) );

saveRDS(lfw_lam_sweep_results, save_filename);
