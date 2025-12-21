library(stringr)
library(tools)

source("./JASA-QuasiLkld_lib.r");


args = commandArgs(TRUE);

if (length(args) < 1) {
  stop(
    paste(
      "No path to data,",
      "file index,",
      "or save directory provided.") );
}

if (length(args) < 3) {
  stop("Missing argument; three arguments expected.");
}

if (length(args) > 3) {
  stop(
    sprintf(
      paste0(
        "run_quasi_for_DEAM_analysis.R takes exactly 3 arguments, ",
        "but %d arguments were provided."),
      length(args) ));
}

DEAM_data_dir = args[1];
file_ix = as.integer(args[2]);
save_dir = args[3];

if (!dir.exists(DEAM_data_dir) ) {
  stop(
    sprintf(
      "The provided data directory %s was not found.",
      DEAM_data_dir) );
}

DEAM_data_filenames = list.files(DEAM_data_dir);
n_DEAM_data_files = length(DEAM_data_filenames);

if ((file_ix < 1) || (file_ix > n_DEAM_data_files) ) {
  stop(
    sprintf(
      "The provided file index %d is invalid. Valid arguments are integers between 1 and %d.",
      file_ix,
      n_DEAM_data_files) );
}

DEAM_data_filename_fx = DEAM_data_filenames[file_ix];
DEAM_data_full_path_fx = file.path(DEAM_data_dir, DEAM_data_filename_fx);

if (file_ext(DEAM_data_filename_fx) != "RData") {
  if (str_length(file_ext(DEAM_data_filename_fx) ) > 0) {
    stop(
      sprintf(
        paste(
          "File of type RData required, ",
          "but filename with extension %s provided.",
          sep=""),
        file_ext(DEAM_data_filename_fx) ));

  } else {
    stop("File of type RData required.");
  }
}

if (!file.exists(save_dir) ) {
  dir.create(save_dir, recursive=TRUE);
}

stopifnot(dir.exists(save_dir) );

DEAM_data = readRDS(DEAM_data_full_path_fx);

stopifnot(class(DEAM_data) == "list");

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
      names(DEAM_data) )) == 0) &&
  (length(
    setdiff(
      names(DEAM_data),
      expected_fields) ) == 0) );

y_train = DEAM_data$y_train;

t_dim = 3;
m_dim_train = length(y_train);

stopifnot(length(DEAM_data$group_ids_train) == m_dim_train);
stopifnot(length(dim(DEAM_data$X_train) ) == t_dim);
stopifnot(dim(DEAM_data$X_train)[t_dim] == m_dim_train);

# Create the fixed effects design matrix X for the training data.
X_train = t(matrix(DEAM_data$X_train, ncol=m_dim_train) );

# The random effects are just random intercepts.
Z_train = matrix(1, m_dim_train, 1);

# Set the group ID's for each observation of the training data.
group_ids_train = DEAM_data$group_ids_train;

set.seed(13);

# Use cross validation to pick an 'a' scalar value from a grid.
a_scal_range = c(1e-4, 10);
n_a_vals = 10;
a_vals = exp(seq(log(a_scal_range[1]), log(a_scal_range[2]), length.out=n_a_vals) );

quasi_cv = CV_quasi_fixed_effects(X_train, y_train, Z_train, group_ids_train, a_vals);
cv_test_mses = rowMeans(as.matrix(quasi_cv[,-1]) );
a_best = max(quasi_cv$a[cv_test_mses == min(cv_test_mses)]);

# Once the best 'a' is chosen, fit the quasi-likelihood model on the training data.
quasi_mod = Fixed_effects_estimation(
  X_train,
  y_train,
  Z_train,
  group_ids_train,
  a_best);

# Estimate the random effects covariance parameter.
resids = as.numeric(y_train - (X_train %*% quasi_mod$beta.hat) );
G_list = list();
G_list[[1]] = diag(x=rep(1, ncol(Z_train) ), ncol(Z_train) );

quasi_rand_eff_var = Varcomp.est(
  resids,
  Z_train,
  G_list,
  a_best,
  group_ids_train);

# Use the estimated random effects covariance parameter to compute the BLUP.
sigma_Z = G_list[[1]] * quasi_rand_eff_var$eta.hat;
blup = compute_blup_quasi(resids, Z_train, sigma_Z, group_ids_train);

y_test = DEAM_data$y_test;

m_dim_test = length(y_test);

stopifnot(length(DEAM_data$group_ids_test) == m_dim_test);
stopifnot(length(dim(DEAM_data$X_test) ) == t_dim);
stopifnot(dim(DEAM_data$X_test)[t_dim] == m_dim_test);

# Create the fixed effects design matrix X for the test data.
X_test = t(matrix(DEAM_data$X_test, ncol=m_dim_test) );

# The random effects are just random intercepts.
Z_test = matrix(1, m_dim_test, 1);

# Set the group ID's for each observation of the test data.
group_ids_test = DEAM_data$group_ids_test;

# Get the quasi-likelihood predictions for the test data using the previously estimated
# beta_hat and the BLUP's.
y_hat_test = predict_quasi(
  X_test,
  Z_test,
  quasi_mod$beta.hat,
  blup,
  group_ids_test);

mspe = mean((y_test - y_hat_test)^2);
pred_r_sqr = 1 - (mspe / mean((y_test - mean(y_test) )^2) );

quasi_info = list(
  quasi_mod=quasi_mod,
  y_test=y_test,
  preds=y_hat_test,
  pred_errs=y_test - y_hat_test,
  mspe=mspe,
  pred_r_sqr=pred_r_sqr);

save_filename = file.path(
  save_dir,
  sprintf(
    "%s_arousal_quasi_results.RData",
    file_path_sans_ext(DEAM_data_filename_fx) ));

saveRDS(quasi_info, save_filename);
