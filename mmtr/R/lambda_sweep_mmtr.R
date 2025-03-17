#' This function performs a mixed model trace regression (MMTR) lambda sweep.
#'
#' @param y Numeric M-vector of responses across all N groups, where M = M1 + M2 + ... + MN.
#' @param X (D+1)-mode fixed effects tensor of dimension P1xP2x...xPDxM, i.e. there is a D-mode
#'   P1xP2x...xPD fixed effects covariate tensor associated with each of the M responses.
#' @param Z (D+1)-mode random effects tensor of dimension Q1xQ2x...xQDxM, i.e. there is a D-mode
#'   Q1xQ2x...xQD random effects covariate tensor associated with each of the M responses.
#' @param group_ids Integer or factor M-vector indicating which group each observation belongs to.
#' @param s_dims_max Integer D-vector indicating the upper bound of the estimated rank of each mode's
#'   random effects covariance matrix.
#' @param n_reps Maximum number of iterations to run the mixed model trace regression AECM algorithm
#'   for each combination of lambda penalty parameters. Defaults to 1000.
#' @param B_lambda_range Minimum and maximum lambda penalty parameters indicating the range of
#'   lambda values over which a parameter sweep will be performed in the estimation of the fixed
#'   effects mean structure B. Defaults to (1e-4, 1).
#' @param L_lambda_range Minimum and maximum lambda penalty parameters indicating the range of
#'   lambda values over which a parameter sweep will be performed in the estimation of the random
#'   effects covariance matrices. Defaults to (1e-4, 1).
#' @param n_B_lambdas Number of different fixed effects lambda penalty parameters to use in lambda
#'   sweep. Defaults to ten.
#' @param n_L_lambdas Number of different random effects lambda penalty parameters to use in lambda
#'   sweep. Defaults to ten.
#' @param comp_ebic Logical value indicating whether or not the EBIC should be computed after the
#'   MMTR model is fit for each combination of lambda penalty parameters. Defaults to TRUE.
#' @param conv_thrsh Convergence threshold used for early stopping of the mixed model trace regression
#'   AECM algorithm based on the change in the objective function from one iteration to the next for
#'   each combination of lambda penalty parameters. Defaults to 1e-3.
#' @param n_proc: Number of processors to use for parallel computation. Defaults to one.
#' @export
lambda_sweep_mmtr = function(
  y, X, Z,
  group_ids,
  s_dims_max,
  n_reps=1000,
  B_lambda_range=c(1e-4, 0.05),
  L_lambda_range=c(1e-4, 1),
  n_B_lambdas=10,
  n_L_lambdas=10,
  comp_ebic=TRUE,
  conv_thrsh=1e-3,
  n_proc=1)
{
  #   Description:
  #
  #     Arguments:
  #                    y: Numeric M-vector of responses across all N groups, where
  #                       M = M1 + M2 + ... + MN.
  #
  #                    X: (D+1)-mode fixed effects tensor of dimension P1xP2x...xPDxM, i.e.
  #                       there is a D-mode P1xP2x...xPD fixed effects covariate tensor
  #                       associated with each of the M responses.
  #
  #                    Z: (D+1)-mode random effects tensor of dimension Q1xQ2x...xQDxM, i.e.
  #                       there is a D-mode Q1xQ2x...xQD random effects covariate tensor
  #                       associated with each of the M responses.
  #
  #            group_ids: Integer or factor M-vector indicating which group each observation
  #                       belongs to.
  #
  #           s_dims_max: Integer D-vector indicating the upper bound of the estimated rank of
  #                       each mode's random effects covariance matrix.
  #
  #               n_reps: Maximum number of iterations to run the mixed model trace regression
  #                       AECM algorithm for each combination of lambda penalty parameters.
  #
  #       B_lambda_range: Minimum and maximum lambda penalty parameters indicating the range of
  #                       lambda values over which a parameter sweep will be performed in the
  #                       estimation of the fixed effects mean structure B.
  #
  #       L_lambda_range: Minimum and maximum lambda penalty parameters indicating the range of
  #                       lambda values over which a parameter sweep will be performed in the
  #                       estimation of the random effects covariance matrices.
  #
  #          n_B_lambdas: Number of different fixed effects lambda penalty parameters to use in
  #                       lambda sweep.
  #
  #          n_L_lambdas: Number of different random effects lambda penalty parameters to use in
  #                       lambda sweep.
  #
  #            comp_ebic: Logical value indicating whether or not the EBIC should be computed after
  #                       the MMTR model is fit for each combination of lambda penalty parameters.
  #
  #           conv_thrsh: Convergence threshold used for early stopping of the mixed model trace
  #                       regression AECM algorithm based on the change in the objective function
  #                       from one iteration to the next for each combination of lambda penalty
  #                       parameters.
  #
  #               n_proc: Number of processors to use for parallel computation.
  #
  #  Return value:

  stopifnot((typeof(y) == "double") || (typeof(y) == "integer") );
  stopifnot(typeof(group_ids) == "integer");
  stopifnot(
    (typeof(s_dims_max) == "double") ||
    (typeof(s_dims_max) == "integer") );

  y = as.vector(y);
  group_ids = as.factor(as.vector(group_ids) );
  s_dims_max = as.vector(s_dims_max);

  if (typeof(y) == "integer") {
    y = as.numeric(y);
  }

  if (typeof(s_dims_max) == "double") {
    s_dims_max = as.integer(s_dims_max);
  }

  stopifnot(all(s_dims_max > 0) );

  stopifnot(is.element("array", class(X) ));
  stopifnot(is.element("array", class(Z) ));

  stopifnot(
    (typeof(B_lambda_range) == "double") ||
    (typeof(B_lambda_range) == "integer") );

  stopifnot(
    (typeof(L_lambda_range) == "double") ||
    (typeof(L_lambda_range) == "integer") );

  stopifnot(length(B_lambda_range) == 2);
  stopifnot(length(L_lambda_range) == 2);

  if (typeof(B_lambda_range) == "integer") {
    B_lambda_range = as.numeric(B_lambda_range);
  }

  if (typeof(L_lambda_range) == "integer") {
    L_lambda_range = as.numeric(L_lambda_range);
  }

  m_dim = length(y);
  n_dim = length(unique(group_ids) );
  t_dim = length(dim(X) ) - 1;

  stopifnot(t_dim == 2);
  stopifnot(length(group_ids) == m_dim);
  stopifnot(length(s_dims_max) == t_dim);
  stopifnot(length(dim(Z) ) == (t_dim + 1) );
  stopifnot(dim(X)[t_dim+1] == m_dim);
  stopifnot(dim(Z)[t_dim+1] == m_dim);

  q_dims = dim(Z)[1:t_dim];
  s_dims = sapply(1:t_dim, function(sx) min(q_dims[sx], s_dims_max[sx]) );

  B_lambdas = exp(
    seq(
      log(B_lambda_range[1]),
      log(B_lambda_range[2]),
      length.out=n_B_lambdas) );

  L_lambdas = exp(
    seq(
      log(L_lambda_range[1]),
      log(L_lambda_range[2]),
      length.out=n_L_lambdas) );

  L_list_init = lapply(
    1:t_dim,
    function(kx) {
      init_mat_kx = matrix(rnorm(q_dims[kx] * s_dims[kx]), q_dims[kx]);
      init_mat_kx_svd = svd(init_mat_kx);
      tcrossprod(init_mat_kx_svd$u, init_mat_kx_svd$v)});

  var_export_vec = c(
    "y",
    "X",
    "Z",
    "group_ids",
    "L_list_init",
    "L_lambdas",
    "comp_ebic",
    "n_reps",
    "conv_thrsh");

  cluster_type = "FORK";

  if (.Platform$OS.type != "unix") {
    cluster_type = "PSOCK";
  }

  n_tries = 10; # Number of tries to call makeCluster().
  wait_time = 5; # In seconds.
  cl = NULL;

  for (ix_try in 1:n_tries) {
    cl = tryCatch(
      parallel::makeCluster(n_proc, type=cluster_type, outfile=""),
      error=\(err) {
        if (ix_try < n_tries) {
          Sys.sleep(wait_time);
        } else {
          stop(err);
        }
      });

    if (!is.null(cl) ) {
      break;
    }
  }

  parallel::clusterEvalQ(cl, library(foreach) );
  parallel::clusterExport(cl, var_export_vec, envir=environment() );
  doParallel::registerDoParallel(cl);

  mmtr_estimates = foreach::foreach(B_lam=B_lambdas) %dopar% {
    n_L_lambdas = length(L_lambdas);
    mmtr_estimates_bx = list();

    for (lx in 1:n_L_lambdas) {
      cat(
        sprintf(
          "Estimating MMTR model with B lambda = %f and L lambda = %f ...\n",
          B_lam,
          L_lambdas[lx]) );

      mmtr_estimates_bx[[lx]] = mmtr::mmtr(
        y, X, Z,
        group_ids,
        L_list_init,
        c(B_lam, L_lambdas[lx]),
        comp_ebic,
        n_reps,
        conv_thrsh,
        1);
    }

    mmtr_estimates_bx
  };

  parallel::stopCluster(cl);

  lambda_sweep_results = list(
    mmtr_estimates=mmtr_estimates,
    B_lambdas=B_lambdas,
    L_lambdas=L_lambdas);

  if (comp_ebic) {
    bic_grid = matrix(NA, n_B_lambdas, n_L_lambdas);

    for (bx in 1:n_B_lambdas) {
      for (lx in 1:n_L_lambdas) {
        bic_grid[bx,lx] = mmtr_estimates[[bx]][[lx]]$ebic;
      }
    }

    bic_min = min(bic_grid);
    ix_bic_min = which(bic_grid == bic_min);
    bx_bic_min = ix_bic_min %% n_B_lambdas;
    bx_bic_min[bx_bic_min == 0] = n_B_lambdas;
    lx_bic_min = floor((ix_bic_min - 1) / n_B_lambdas) + 1;

    ix_bic_min = which(bx_bic_min == max(bx_bic_min) );
    ix_bic_min = which(lx_bic_min[ix_bic_min] == max(lx_bic_min[ix_bic_min]) );
    bx_bic_min = bx_bic_min[ix_bic_min];
    lx_bic_min = lx_bic_min[ix_bic_min];

    lambda_sweep_results$mmtr_mod_min_ebic = mmtr_estimates[[bx_bic_min]][[lx_bic_min]];
  }

  lambda_sweep_results
};
