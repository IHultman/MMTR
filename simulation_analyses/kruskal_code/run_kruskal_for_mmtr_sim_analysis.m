function [] = kruskal_simulation_analysis( ...
  sim_data_dir, ...
  results_save_dir, ...
  ix_sim_set, ...
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

  sim_data_dir = char(sim_data_dir);
  results_save_dir = char(results_save_dir);

  assert( ...
    isfolder(sim_data_dir), ...
    ['The provided directory ', sim_data_dir, ' was not found.']);

  if ~isfolder(results_save_dir)
    mkdir(results_save_dir);

    assert(isfolder(results_save_dir) );
  end

  if isstring(ix_sim_set) || ischar(ix_sim_set)
    ix_sim_set = str2num(ix_sim_set);
  end

  if isstring(rank) || ischar(rank)
    rank = str2num(rank);
  end

  sim_sets = split(ls(sim_data_dir) );
  sim_sets = sort(sim_sets(1:(end-1) ));
  n_sim_sets = numel(sim_sets);

  assert( ...
    (ix_sim_set >= 1) && (ix_sim_set <= n_sim_sets), ...
    sprintf( ...
      'The provided file index %d is invalid. Valid arguments are integers between 1 and %d.', ...
      ix_sim_set, ...
      n_sim_sets) );

  sim_data_dir_sx = [sim_data_dir, '/', sim_sets{ix_sim_set}, '/'];

  assert( ...
    isfolder(sim_data_dir_sx), ...
    ['The file ', sim_data_dir_sx, ' is not a directory.']);

  sim_filenames_sx = split(ls(sim_data_dir_sx) );
  sim_filenames_sx = sim_filenames_sx(1:(end-1) );
  n_sim_filenames = numel(sim_filenames_sx);

  disp(sim_filenames_sx);

  save_dir_sx = [results_save_dir, '/', sim_sets{ix_sim_set}, '/'];

  if ~isfolder(save_dir_sx)
    mkdir(save_dir_sx);

    assert(isfolder(save_dir_sx) );
  end

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

  dist = 'normal';
  pentype = 'enet';
  penparam = 1;

  disp(['Fitting Kruskal model for simulation ', sim_sets{ix_sim_set}, ' ...']);

  for fx = 1:n_sim_filenames
    disp( ...
      ['Fitting Kruskal model for simulation scenario ', ...
       sim_filenames_sx{fx}, '  ', num2str(fx), '/', num2str(n_sim_filenames), ' ...']);

    sim_data_full_path_fx = [sim_data_dir_sx, sim_filenames_sx{fx}];
    sim_data = load(sim_data_full_path_fx);

    t_dim = 2;
    m_dim = size(sim_data.X, t_dim + 1);
    p_dims = size(sim_data.X, 1:t_dim);
    rank = min(rank, min(p_dims) );
    X_int = ones(m_dim, 1);

    assert(length(sim_data.y) == m_dim);
  
    err_mse = kruskal_cv( ...
      X_int, ...
      sim_data.X, ...
      sim_data.y, ...
      rank, ...
      dist, ...
      lambdas, ...
      pentype, penparam, ...
      false);

    ix_min_mse = max(find(err_mse == min(err_mse) ));
    best_lambda = lambdas(ix_min_mse);

    [~, beta_t_init] = kruskal_reg(X_int, sim_data.X, sim_data.y, rank, dist);

    [beta_0, beta_t, ~] = kruskal_sparsereg( ...
      X_int, ...
      sim_data.X, ...
      sim_data.y, ...
      rank, ...
      dist, ...
      best_lambda, ...
      pentype, penparam, ...
      'B0', beta_t_init);

    y_fitted = ( ...
      (X_int * beta_0) + ...
      double(ttt(tensor(beta_t), tensor(sim_data.X), 1:2) ));

    results_struct = struct();
    results_struct.b0_hat = beta_0;
    results_struct.B_hat = double(beta_t);
    results_struct.y = sim_data.y;
    results_struct.y_fitted = y_fitted;
    results_struct.resids = sim_data.y - y_fitted;
    results_struct.tau2_hat = mse(y_fitted, sim_data.y);

    [~, sim_filename_no_ext_fx] = fileparts(sim_filenames_sx{fx});
    save_filename = [save_dir_sx, sim_filename_no_ext_fx, '_cv_results.mat'];

    save(save_filename, '-struct', 'results_struct');
  end
end
