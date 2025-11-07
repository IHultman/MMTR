library(dplyr)
library(R.matlab)


ascertain_data_fname = "./ASCERTAIN_dataset.mat";

r_save_dir = "./mmtr_ASCERTAIN_train_test_sets";

if (!file.exists(r_save_dir) ) {
  dir.create(r_save_dir);
}

stopifnot(dir.exists(r_save_dir) );

matlab_save_dir = "./kruskal_ASCERTAIN_train_test_sets";

if (!file.exists(matlab_save_dir) ) {
  dir.create(matlab_save_dir);
}

stopifnot(dir.exists(matlab_save_dir) );

ascertain_data = readMat(ascertain_data_fname);

n_features = 11;
n_chans = 8;
n_video = 36;
n_subj = 58;

n_tot_obs = n_video * n_subj;

X = ascertain_data$X;
group_ids = as.vector(ascertain_data$group.ids);
y = ascertain_data$y;

stopifnot(all(dim(X) == c(n_features, n_chans, n_tot_obs) ));
stopifnot(all(dim(y) == c(n_tot_obs, 5) ));

group_ids_check = table(group_ids);

stopifnot(all(group_ids_check == n_video) );
stopifnot(length(group_ids_check) == n_subj);
stopifnot(all(sort(group_ids) == group_ids) );

group_ixs_ranges = t(
  data.frame(group_id=group_ids, ix=1:length(group_ids) ) %>%
  group_by(group_id) %>%
  summarize(ix_start=min(ix), ix_end=max(ix), .groups="drop") %>%
  select(ix_start, ix_end) %>%
  as.matrix() );

set.seed(13);

n_test_eeg_per_grp = 5;
n_train_test_sets = 10;

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

  y_train = y[train_ixs,];
  rownames(y_train) = NULL;

  y_test = y[test_ixs,];
  rownames(y_test) = NULL;

  group_ids_train = as.factor(group_ids[train_ixs]);
  group_ids_test = as.factor(group_ids[test_ixs]);

  stopifnot(dim(X_train)[3] == nrow(y_train) );
  stopifnot(length(group_ids_train) == nrow(y_train) );
  stopifnot(dim(X_test)[3] == nrow(y_test) );
  stopifnot(length(group_ids_test) == nrow(y_test) );

  save_filename = paste0("ASCERTAIN_train_test_set_", sx);

  saveRDS(
    list(
      X_train=X_train,
      X_test=X_test,
      y_train=y_train,
      y_test=y_test,
      group_ids_train=group_ids_train,
      group_ids_test=group_ids_test),
    file.path(r_save_dir, paste0(save_filename, ".RData") ));

  writeMat(
    file.path(matlab_save_dir, paste0(save_filename, ".mat") ),
    X_train=X_train,
    X_test=X_test,
    y_train=as.matrix(y_train),
    y_test=as.matrix(y_test),
    group_ids_train=as.character(group_ids_train),
    group_ids_test=as.character(group_ids_test) );
}
