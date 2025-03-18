function [] = kruskal_lfw_analysis( ...
  lfw_data_dir, ...
  results_save_dir, ...
  file_ix, ...
  attribute_of_interest, ...
  rank)

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

  if isstring(file_ix) || ischar(file_ix)
    file_ix = str2num(file_ix);
  end

  if isstring(rank) || ischar(rank)
    rank = str2num(rank);
  end

  lfw_data_filenames = split(ls(lfw_data_dir) );
  lfw_data_filenames = sort(lfw_data_filenames(1:(end-1) ));
  n_lfw_data_files = numel(lfw_data_filenames);

  assert( ...
    (file_ix >= 1) && (file_ix <= n_lfw_data_files), ...
    sprintf( ...
      'The provided file index %d is invalid. Valid arguments are integers between 1 and %d.', ...
      file_ix, ...
      n_lfw_data_files) );

  lfw_data_filename_fx = lfw_data_filenames{file_ix};
  lfw_data_full_path_fx = [lfw_data_dir, '/', lfw_data_filename_fx];
  lfw_data = load(lfw_data_full_path_fx);

  lfw_attributes_all = fieldnames(lfw_data.y_train);
  ix_match = strcmp(attribute_of_interest, lfw_attributes_all);

  assert(sum(ix_match) > 0, ['Attribute ', attribute_of_interest, ' not found.']);
  assert( ...
    sum(ix_match) == 1, ...
    ['Multiple attribute variables match ', attribute_of_interest, '.']);

  disp(['Fitting Kruskal model for attribute ', attribute_of_interest ' ...']);

  M_train = lfw_data.X_train;
  M_test = lfw_data.X_test;

  y_train = lfw_data.y_train.(attribute_of_interest);
  y_test = lfw_data.y_test.(attribute_of_interest);
  
  t_dim = length(size(M_train) );
  p_dims = size(M_train, 1:(t_dim-1) );
  m_dim_train = size(M_train, t_dim);
  m_dim_test = size(M_test, t_dim);
  
  assert(length(y_train) == m_dim_train);
  assert(length(y_test) == m_dim_test);
  
  X_int_train = ones(m_dim_train, 1);
  X_int_test = ones(m_dim_test, 1);

  lambda_range = [1e-2, 10];
  n_lambdas = 10;

  if (n_lambdas > 1)
    log_lambda_range = log(lambda_range);
    log_lam_step = abs(diff(log_lambda_range) ) / (n_lambdas - 1);
    log_lambdas = log_lambda_range(1):log_lam_step:log_lambda_range(2);
    lambdas = exp(log_lambdas);
  else
    lambdas = lambda_range(1);
  end

  rank = min(rank, min(p_dims) );
  dist = 'normal';
  pentype = 'enet';
  penparam = 1;

  disp(sprintf("Performing cross validation with rank %d mean structure ...", rank) );

  err_mse = kruskal_cv( ...
    X_int_train, ...
    M_train, ...
    y_train, ...
    rank, ...
    dist, ...
    lambdas, ...
    pentype, penparam, ...
    false);

  min_mse_ix = max(find(err_mse == min(err_mse) ));
  best_lambda = lambdas(min_mse_ix);

  [~, beta_t_init] = kruskal_reg(X_int_train, M_train, y_train, rank, dist);

  [beta_0, beta_t, ~] = kruskal_sparsereg( ...
    X_int_train, ...
    M_train, ...
    y_train, ...
    rank, ...
    dist, ...
    best_lambda, ...
    pentype, penparam, ...
    'B0', beta_t_init);

  y_test_preds = ( ...
    (X_int_test * beta_0) + ...
    double(ttt(tensor(beta_t), tensor(M_test), 1:2) ));

  results_struct = struct();
  results_struct.b0_hat = beta_0;
  results_struct.B_hat = double(tensor(beta_t) );
  results_struct.preds = y_test_preds;
  results_struct.pred_errs = y_test - y_test_preds;
  results_struct.mspe = mean((y_test - y_test_preds).^2);
  results_struct.pred_r_sqr = 1 - (results_struct.mspe / mean((y_test - mean(y_test) ).^2) );

  [~, lfw_filename_no_ext_fx] = fileparts(lfw_data_filename_fx);
  save_filename = [ ...
    results_save_dir, ...
    sprintf( ...
      '/%s_%s_rank_%d_kruskal_cv_results.mat', ...
      lfw_filename_no_ext_fx, ...
      attribute_of_interest, ...
      rank)];

  save(save_filename, '-struct', 'results_struct');
end
