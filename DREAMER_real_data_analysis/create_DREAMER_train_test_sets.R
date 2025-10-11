library(R.matlab)


processed_dreamer_filename = "./processed_DREAMER_by_freq_bands.mat";

mmtr_save_dir = "./mmtr_DREAMER_train_test_sets";

if (!file.exists(mmtr_save_dir) ) {
  dir.create(mmtr_save_dir);
}

stopifnot(dir.exists(mmtr_save_dir) );

kruskal_save_dir = "./kruskal_DREAMER_train_test_sets";

if (!file.exists(kruskal_save_dir) ) {
  dir.create(kruskal_save_dir);
}

stopifnot(dir.exists(kruskal_save_dir) );

dreamer_data = readMat(processed_dreamer_filename);

p_dims = dim(dreamer_data$eeg.stim.logpower)[1:2];
m_dim = dim(dreamer_data$eeg.stim.logpower)[3];
n_dim = dim(dreamer_data$eeg.stim.logpower)[4];
n_tot_obs = m_dim * n_dim;

X = array(dreamer_data$eeg.stim.logpower, dim=c(p_dims, n_tot_obs) );

score_labels = unlist(dreamer_data$score.labels);

stopifnot(dim(dreamer_data$scores)[1] == m_dim);
stopifnot(dim(dreamer_data$scores)[2] == length(score_labels) );
stopifnot(dim(dreamer_data$scores)[3] == n_dim);

y_df = as.data.frame(matrix(aperm(dreamer_data$scores, c(1, 3, 2) ), nrow=n_tot_obs) );
colnames(y_df) = score_labels;

set.seed(13);

n_test_eeg_per_grp = 4;
n_train_test_sets = 10;

group_ids = rep(1:n_dim, each=m_dim);

group_ixs_ranges = rbind(
  seq(1, n_tot_obs, m_dim),
  seq(m_dim, n_tot_obs, m_dim) );

for (sx in 1:n_train_test_sets) {
  test_ixs = sort(
    as.numeric(
      apply(
        group_ixs_ranges,
        MARGIN=2,
        \(col_ix) sample(col_ix[1]:col_ix[2], n_test_eeg_per_grp) )));

  train_ixs = setdiff(1:n_tot_obs, test_ixs);

  X_train = X[,,train_ixs];
  X_test = X[,,test_ixs];

  y_train = y_df[train_ixs,];
  rownames(y_train) = NULL;

  y_test = y_df[test_ixs,];
  rownames(y_test) = NULL;

  group_ids_train = as.factor(group_ids[train_ixs]);
  group_ids_test = as.factor(group_ids[test_ixs]);

  stopifnot(dim(X_train)[3] == nrow(y_train) );
  stopifnot(length(group_ids_train) == nrow(y_train) );
  stopifnot(dim(X_test)[3] == nrow(y_test) );
  stopifnot(length(group_ids_test) == nrow(y_test) );

  save_filename = paste0("DREAMER_train_test_set_", sx);

  saveRDS(
    list(
      X_train=X_train,
      X_test=X_test,
      y_train=y_train,
      y_test=y_test,
      group_ids_train=group_ids_train,
      group_ids_test=group_ids_test),
    file.path(mmtr_save_dir, paste0(save_filename, ".RData") ));

  writeMat(
    file.path(kruskal_save_dir, paste0(save_filename, ".mat") ),
    X_train=X_train,
    X_test=X_test,
    y_train=as.matrix(y_train),
    y_test=as.matrix(y_test),
    score_labels=score_labels,
    group_ids_train=as.character(group_ids_train),
    group_ids_test=as.character(group_ids_test) );
}
