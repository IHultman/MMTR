function [err_mse, err_se, auc_mse, auc_se] = kruskal_cv( ...
  X, M, y, ...
  rnk, ...
  dist, ...
  lambdas, ...
  pentype, penparam, ...
  include_AUC)
%{
DESCRIPTION:
  Get the 10-fold cross validation MSE's (optionally with AUC's) for the rank
  $rnk estimation of the parameter tensor corresponding to the input tensor
  observations $M at input lambda(s).

PARAMETERS:
            X:
            M:
            y:
          rnk:
         dist:
      lambdas:
      pentype:
     penparam:
  include_AUC:

RETURN VALUE:
  10-fold cross validation MSE's (optionally with AUC's) corresponding to the
  provided lambda value(s).
%}

  assert(length(size(X) ) == 2);
  assert(length(size(M) ) > 2);
  assert(size(X, 1) == length(y) );
  assert(size(M, length(size(M) )) == length(y) );

  X = double(X);
  M = double(M);
  y = double(y);

  n_lams = length(lambdas);

  n_cv_grps = 10;
  n_tot = length(y);
  ns_vec = repelem(floor(n_tot / n_cv_grps), n_cv_grps);
  n_rmnng = n_tot - sum(ns_vec);

  if (n_rmnng > 0)
    ns_vec(1:n_rmnng) = ns_vec(1:n_rmnng) + 1;
  end

  rndm_ixs = randsample(n_tot, n_tot);
  grp_ixs = [ ...
    cumsum([1, ns_vec(1:(end-1) )])', ...
    cumsum(ns_vec)'];

  errs = repelem(nan, n_lams, n_cv_grps);
  aucs = repelem(nan, n_lams, n_cv_grps);

  for gx = 1:n_cv_grps
    tst_ixs = sort(rndm_ixs(grp_ixs(gx,1):grp_ixs(gx,2) ));
    trn_ixs = setdiff(1:n_tot, tst_ixs)';

    [~, beta_t_init] = kruskal_reg( ...
      X(trn_ixs,:), ...
      M(:,:,trn_ixs), ...
      y(trn_ixs,1), ...
      rnk, ...
      dist);

    for lx = 1:n_lams
      [beta_0, beta_t, ~] = kruskal_sparsereg( ...
        X(trn_ixs,:), ...
        M(:,:,trn_ixs), ...
        y(trn_ixs,1), ...
        rnk, ...
        dist, ...
        lambdas(lx), ...
        pentype, penparam, ...
        'B0', beta_t_init);

      y_preds = ( ...
        (X(tst_ixs,:) * beta_0) + ...
        double( ...
          ttt( ...
            tensor(beta_t), ...
            tensor(M(:,:,tst_ixs) ), ...
            1:2) ));

      if (dist == "binomial")
        y_preds = 1 ./ (1 + exp(-y_preds) );

        if (include_AUC)
          [~, ~, ~, auc_val] = perfcurve(y(tst_ixs,1), y_preds, 1);
          aucs(lx,gx) = auc_val;
        end
      end

      errs(lx,gx) = mse(y_preds, y(tst_ixs,1) );
    end
  end

  err_mse = mean(errs, 2);
  err_se = std(errs')' / sqrt(n_cv_grps);
  
  auc_mse = mean(aucs, 2);
  auc_se = std(aucs')' / sqrt(n_cv_grps);
end

