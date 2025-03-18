library(lme4)
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

p_dims = dim(lfw_data$X_train)[1:2];
P_dim = prod(p_dims);

lfw_train_df = data.frame(
  y=lfw_data$y_train[[lfw_attribute]],
  group=lfw_data$group_labels_train,
  t(matrix(lfw_data$X_train, P_dim) ));

lfw_test_df = data.frame(
  y=lfw_data$y_test[[lfw_attribute]],
  group=lfw_data$group_labels_test,
  t(matrix(lfw_data$X_test, P_dim) ));

lfw_lmer_formula = as.formula(
  sprintf(
    "y ~ %s + (1 | group)",
    paste0(
      setdiff(colnames(lfw_train_df), c("y", "group") ),
      collapse=" + ") ));

set.seed(13);

lfw_lmer_mod = lmer(lfw_lmer_formula, data=lfw_train_df);
lfw_test_preds = as.numeric(predict(lfw_lmer_mod, lfw_test_df) );

if (!dir.exists(save_dir) ) {
  dir.create(save_dir, recursive=TRUE);

  stopifnot(dir.exists(save_dir) );
}

save_filename = paste0(
  save_dir,
  sprintf(
    "/%s_%s_lme4_results.RData",
    file_path_sans_ext(lfw_data_filename_fx),
    lfw_attribute) );

saveRDS(
  list(
    lmer_mod=lfw_lmer_mod,
    preds=lfw_test_preds,
    pred_errs=lfw_test_df$y - lfw_test_preds,
    mspe=mean((lfw_test_df$y - lfw_test_preds)^2) ),
  save_filename);
