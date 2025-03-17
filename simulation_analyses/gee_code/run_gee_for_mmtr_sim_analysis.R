function [] = gee_simulation_analysis( ...
  sim_data_dir, ...
  results_save_dir, ...
  ix_sim_set)

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
  %covar_type = 'unstructured';
  covar_type = 'equicorr';
  pentype = 'enet';
  penparam = 1;

  disp(['Fitting Tensor GEE model for simulation ', sim_sets{ix_sim_set}, ' ...']);

  for fx = 1:n_sim_filenames
    disp( ...
      ['Fitting Tensor GEE model for simulation scenario ', ...
       sim_filenames_sx{fx}, '  ', num2str(fx), '/', num2str(n_sim_filenames), ' ...']);

    sim_data_full_path_fx = [sim_data_dir_sx, sim_filenames_sx{fx}];
    sim_data = load(sim_data_full_path_fx);

    vec_n_obs = groupcounts(sim_data.group_ids);
    time = arrayfun(@(n_obs) 1:n_obs, vec_n_obs, 'UniformOutput', false);
    time = [time{:}];

    if isrow(time)
      time = time';
    end

    t_dim = 2;
    m_dim = size(sim_data.X, t_dim + 1);
    P_dim = prod(size(sim_data.X, 1:t_dim) );
    X = reshape(sim_data.X, P_dim, m_dim)';

    err_mse = gee_cv( ...
      X, sim_data.y, ...
      sim_data.group_ids, ...
      time, ...
      dist, ...
      covar_type, ...
      lambdas, ...
      pentype, ...
      penparam);

    ix_min_mse = max(find(err_mse == min(err_mse) ));
    cv_lambda = lambdas(ix_min_mse);

    [betahat, alphahat, stats] = gee_sparsereg( ...
      findgroups(sim_data.group_ids), ...
      time, ...
      X, ...
      sim_data.y, ...
      dist, ...
      covar_type, ...
      cv_lambda, ...
      'penalty', pentype, ...
      'penparam', penparam);

    y_fitted = X * betahat;

    if strcmp(covar_type, 'equicorr')
      m_gx_vec = groupcounts(sim_data.group_ids);
      m_gx = m_gx_vec(1);

      assert(all(m_gx_vec == m_gx) );

      alphahat = alphahat * ones(m_gx, m_gx);
      ixs_diag = 1:(m_gx+1):(m_gx^2);
      alphahat(ixs_diag) = 1;
    end

    results_struct = struct();
    results_struct.betahat = betahat;
    results_struct.alphahat = alphahat;
    results_struct.y = sim_data.y;
    results_struct.y_fitted = y_fitted;
    results_struct.resids = sim_data.y - y_fitted;
    results_struct.covar_type = covar_type;

    [~, sim_filename_no_ext_fx] = fileparts(sim_filenames_sx{fx});
    save_filename = [save_dir_sx, sim_filename_no_ext_fx, '_', char(covar_type), '_gee_cv_results.mat'];

    save(save_filename, '-struct', 'results_struct');
  end
end
