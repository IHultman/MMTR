library(dplyr)
library(R.matlab)


DEAM_data_fname = file.path(
  "./DEAM_raw_data",
  "DEAM_preprocessed_data_averaged_across_rater_arousal_response.RData");

r_save_dir = "./mmtr_DEAM_arousal_train_test_sets";

if (!file.exists(r_save_dir) ) {
  dir.create(r_save_dir);
}

stopifnot(dir.exists(r_save_dir) );

matlab_save_dir = "./kruskal_DEAM_arousal_train_test_sets";

if (!file.exists(matlab_save_dir) ) {
  dir.create(matlab_save_dir);
}

stopifnot(dir.exists(matlab_save_dir) );

DEAM_data = readRDS(DEAM_data_fname);

p1_dim = dim(DEAM_data$X)[1];
p2_dim = dim(DEAM_data$X)[2];
m_dim = nrow(DEAM_data$group_info);

stopifnot(dim(DEAM_data$X)[3] == m_dim);

ixs_order = order(DEAM_data$group_info$song_id, DEAM_data$group_info$time_bin);
DEAM_data$X = DEAM_data$X[,,ixs_order];
DEAM_data$group_info = DEAM_data$group_info[ixs_order,];

X = DEAM_data$X;
group_ids = DEAM_data$group_info$song_id;
y = DEAM_data$group_info$y;

group_ixs_ranges = t(
  data.frame(group_id=group_ids, ix=1:length(group_ids) ) %>%
  group_by(group_id) %>%
  summarize(ix_start=min(ix), ix_end=max(ix), .groups="drop") %>%
  select(ix_start, ix_end) %>%
  as.matrix() );

set.seed(13);

n_test_obs_per_grp = 7;
n_train_test_sets = 10;

for (sx in 1:n_train_test_sets) {
  test_ixs = sort(
    as.numeric(
      apply(
        group_ixs_ranges,
        MARGIN=2,
        \(col_ix) sample(col_ix[1]:col_ix[2], n_test_obs_per_grp) )));

  train_ixs = setdiff(1:m_dim, test_ixs);

  X_train = X[,,train_ixs];
  X_test = X[,,test_ixs];

  y_train = y[train_ixs];
  y_test = y[test_ixs];

  group_ids_train = factor(group_ids[train_ixs], levels=unique(group_ids) );
  group_ids_test = factor(group_ids[test_ixs], levels=unique(group_ids) );

  stopifnot(dim(X_train)[3] == length(y_train) );
  stopifnot(length(group_ids_train) == length(y_train) );
  stopifnot(dim(X_test)[3] == length(y_test) );
  stopifnot(length(group_ids_test) == length(y_test) );

  save_filename = paste0("DEAM_train_test_set_", sx);

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
    group_ids_train=as.numeric(group_ids_train),
    group_ids_test=as.numeric(group_ids_test) );
}
