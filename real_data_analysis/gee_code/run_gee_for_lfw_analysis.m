function [] = gee_lfw_analysis( ...
  lfw_data_dir, ...
  results_save_dir, ...
  ix_file, ...
  attribute_of_interest)

  kruskal_lib_dir = [ ...
    '/Users', ...
    '/ikhultman', ...
    '/mmtr_paper_final', ...
    '/mmtr_paper_final', ...
    '/kruskal_analysis', ...
    '/kruskal_code'];

  addpath(genpath([kruskal_lib_dir, '/SparseReg']) );
  addpath(genpath([kruskal_lib_dir, '/TensorReg']) );
  addpath(genpath([kruskal_lib_dir, '/tensor_toolbox']) );

  lfw_data_dir = char(lfw_data_dir);
  results_save_dir = char(results_save_dir);
  attribute_of_interest = char(attribute_of_interest);

  assert( ...
    isfolder(lfw_data_dir), ...
    ['The provided directory ', lfw_data_dir, ' was not found.']);

  if ~isfolder(results_save_dir)
    mkdir(results_save_dir);

    assert(isfolder(results_save_dir) );
  end

  if isstring(ix_file) || ischar(ix_file)
    ix_file = str2num(ix_file);
  end

  lfw_data_filenames = split(ls(lfw_data_dir) );
  lfw_data_filenames = sort(lfw_data_filenames(1:(end-1) ));
  n_lfw_data_files = numel(lfw_data_filenames);

  assert( ...
    (ix_file >= 1) && (ix_file <= n_lfw_data_files), ...
    sprintf( ...
      'The provided file index %d is invalid. Valid arguments are integers between 1 and %d.', ...
      ix_file, ...
      n_lfw_data_files) );

  lfw_data_filename_fx = lfw_data_filenames{ix_file};
  lfw_data_full_path_fx = [lfw_data_dir, '/', lfw_data_filename_fx];
  lfw_data = load(lfw_data_full_path_fx);

  lfw_attributes_all = fieldnames(lfw_data.y_train);
  ix_match = strcmp(attribute_of_interest, lfw_attributes_all);

  assert(sum(ix_match) > 0, ['Attribute ', attribute_of_interest, ' not found.']);
  assert( ...
    sum(ix_match) == 1, ...
    ['Multiple attribute variables match ', attribute_of_interest, '.']);

  disp(['Fitting GEE model for attribute ', attribute_of_interest ' ...']);

  lambda_range = [1e-4, 10];
  n_lambdas = 10;

  if (n_lambdas > 1)
    log_lambda_range = log(lambda_range);
    log_lam_step = abs(diff(log_lambda_range) ) / (n_lambdas - 1);
    log_lambdas = log_lambda_range(1):log_lam_step:log_lambda_range(2);
    lambdas = exp(log_lambdas);
  else
    lambdas = lambda_range(1);
  end

  dist = 'normal';
  covar_type = 'equicorr';
  pentype = 'enet';
  penparam = 1;

  t_dim = 2;
  m_dim_train = size(lfw_data.X_train, t_dim + 1);
  m_dim_test = size(lfw_data.X_test, t_dim + 1);
  P_dim = prod(size(lfw_data.X_train, 1:t_dim) );

  X_train = reshape(lfw_data.X_train, P_dim, m_dim_train)';
  X_test = reshape(lfw_data.X_test, P_dim, m_dim_test)';

  y_train = lfw_data.y_train.(attribute_of_interest);
  y_test = lfw_data.y_test.(attribute_of_interest);

  group_labels_train = lfw_data.group_labels_train;
  group_labels_test = lfw_data.group_labels_test;

  vec_n_train_obs = groupcounts(group_labels_train);
  time_train = arrayfun(@(n_obs) 1:n_obs, vec_n_train_obs, 'UniformOutput', false);
  time_train = [time_train{:}];

  if isrow(time_train)
    time_train = time_train';
  end

  vec_n_test_obs = groupcounts(group_labels_test);
  time_test = arrayfun(@(n_obs) 1:n_obs, vec_n_test_obs, 'UniformOutput', false);
  time_test = [time_test{:}];

  if isrow(time_test)
    time_test = time_test';
  end

  assert(length(y_train) == m_dim_train);
  assert(length(y_test) == m_dim_test);
  assert(length(group_labels_train) == m_dim_train);
  assert(length(group_labels_test) == m_dim_test);
  assert(length(time_train) == m_dim_train);
  assert(length(time_test) == m_dim_test);

  err_mse = gee_cv( ...
    X_train, ...
    y_train, ...
    group_labels_train, ...
    time_train, ...
    dist, ...
    covar_type, ...
    lambdas, ...
    pentype, ...
    penparam);

  ix_min_mse = max(find(err_mse == min(err_mse) ));
  cv_lambda = lambdas(ix_min_mse);

  [betahat, alphahat, stats] = gee_sparsereg( ...
    findgroups(group_labels_train), ...
    time_train, ...
    X_train, ...
    y_train, ...
    dist, ...
    covar_type, ...
    cv_lambda, ...
    'penalty', pentype, ...
    'penparam', penparam);

  resids = y_train - (X_train * betahat);

  y_preds = gee_equicorr_predict( ...
    resids,  ...
    group_labels_train, ...
    betahat, ...
    alphahat, ...
    X_test, ...
    group_labels_test);

  results_struct = struct();
  results_struct.betahat = betahat;
  results_struct.alphahat = alphahat;
  results_struct.preds = y_preds;
  results_struct.pred_errs = y_test - y_preds;
  results_struct.mspe = mean((y_test - y_preds).^2);
  results_struct.pred_r_sqr = 1 - (results_struct.mspe / mean((y_test - mean(y_test) ).^2) );

  [~, lfw_filename_no_ext_fx] = fileparts(lfw_data_filename_fx);
  save_filename = [ ...
    results_save_dir, ...
    sprintf( ...
      '/%s_%s_gee_cv_results.mat', ...
      lfw_filename_no_ext_fx, ...
      attribute_of_interest)];

  save(save_filename, '-struct', 'results_struct');
end
