library(dplyr)
library(stringr)
library(tidyr)


# You can change this path to the valence file if you want valence responses.
annotations_filename = file.path(
  "./DEAM_raw_data",
  "DEAM_Annotations",
  "annotations",
  "annotations averaged per song",
  "dynamic (per second annotations)",
  "valence.csv");

features_dir = "./DEAM_raw_data/features/features";

save_filename = file.path(
  "./DEAM_train_test_sets_in_order_times",
  "DEAM_preprocessed_data_averaged_across_rater_valence_response.RData");

# Use only the songs whose ID's are between 1001 and 2000. The data for this set of songs
# are all of length approximately 43 seconds.
song_ids = 1001:2000;
n_songs = length(song_ids);

annotations = (
  read.csv(annotations_filename, header=TRUE) %>%
  filter(is.element(song_id, song_ids) ));

ixs_cols = which(colSums(!is.na(annotations) ) > 0);
annotations = annotations[,ixs_cols];

time_labels_init = str_match(names(annotations), "^sample_[0-9]+ms$");
time_labels_init = time_labels_init[!is.na(time_labels_init)];

# Create a map from the columns of the annotations table to the time in seconds to which
# they correspond.
time_labels_map_init = data.frame(
  time_label=time_labels_init,
  time=as.numeric(str_extract(time_labels_init, "(?<=^sample_)[0-9]+(?=ms$)") ) / 1000);

feat_grps_patterns = c(
#  "audSpec_Rfilt_sma%s.0._%s",
#  "audSpec_Rfilt_sma%s.1._%s",
#  "audSpec_Rfilt_sma%s.2._%s",
#  "audSpec_Rfilt_sma%s.3._%s",
#  "audSpec_Rfilt_sma%s.4._%s",
#  "audSpec_Rfilt_sma%s.5._%s",
#  "audSpec_Rfilt_sma%s.6._%s",
#  "audSpec_Rfilt_sma%s.7._%s",
#  "audSpec_Rfilt_sma%s.8._%s",
#  "audSpec_Rfilt_sma%s.9._%s",
#  "audSpec_Rfilt_sma%s.10._%s",
#  "audSpec_Rfilt_sma%s.11._%s",
#  "audSpec_Rfilt_sma%s.12._%s",
#  "audSpec_Rfilt_sma%s.13._%s",
#  "audSpec_Rfilt_sma%s.14._%s",
#  "audSpec_Rfilt_sma%s.15._%s",
#  "audSpec_Rfilt_sma%s.16._%s",
#  "audSpec_Rfilt_sma%s.17._%s",
#  "audSpec_Rfilt_sma%s.18._%s",
#  "audSpec_Rfilt_sma%s.19._%s",
#  "audSpec_Rfilt_sma%s.20._%s",
#  "audSpec_Rfilt_sma%s.21._%s",
#  "audSpec_Rfilt_sma%s.22._%s",
#  "audSpec_Rfilt_sma%s.23._%s",
#  "audSpec_Rfilt_sma%s.24._%s",
#  "audSpec_Rfilt_sma%s.25._%s",
  "pcm_fftMag_fband250.650_sma%s_%s",
  "pcm_fftMag_fband1000.4000_sma%s_%s",
  "pcm_fftMag_spectralRollOff25.0_sma%s_%s",
  "pcm_fftMag_spectralRollOff50.0_sma%s_%s",
  "pcm_fftMag_spectralRollOff75.0_sma%s_%s",
  "pcm_fftMag_spectralRollOff90.0_sma%s_%s",
  "pcm_fftMag_spectralFlux_sma%s_%s",
  "pcm_fftMag_spectralCentroid_sma%s_%s",
  "pcm_fftMag_spectralEntropy_sma%s_%s",
  "pcm_fftMag_spectralVariance_sma%s_%s",
  "pcm_fftMag_spectralSkewness_sma%s_%s",
  "pcm_fftMag_spectralKurtosis_sma%s_%s",
  "pcm_fftMag_spectralSlope_sma%s_%s",
  "pcm_fftMag_psySharpness_sma%s_%s",
  "pcm_fftMag_spectralHarmonicity_sma%s_%s",
  "pcm_fftMag_mfcc_sma%s.1._%s",
  "pcm_fftMag_mfcc_sma%s.2._%s",
  "pcm_fftMag_mfcc_sma%s.3._%s",
  "pcm_fftMag_mfcc_sma%s.4._%s",
  "pcm_fftMag_mfcc_sma%s.5._%s",
  "pcm_fftMag_mfcc_sma%s.6._%s",
  "pcm_fftMag_mfcc_sma%s.7._%s",
  "pcm_fftMag_mfcc_sma%s.8._%s",
  "pcm_fftMag_mfcc_sma%s.9._%s",
  "pcm_fftMag_mfcc_sma%s.10._%s",
  "pcm_fftMag_mfcc_sma%s.11._%s",
  "pcm_fftMag_mfcc_sma%s.12._%s",
  "pcm_fftMag_mfcc_sma%s.13._%s",
  "pcm_fftMag_mfcc_sma%s.14._%s");

feat_grps = lapply(
  list(c("", "amean"), c("", "stddev"), c("_de","amean"), c("_de", "stddev") ),
  \(grp_lbl) do.call(sprintf, c(list(feat_grps_patterns), grp_lbl) ));

feat_grp_labels = unlist(feat_grps);

n_feat_grps = length(feat_grps);
len_feat_grps = length(feat_grps_patterns);

combined_datasets = list();

for (sx in 1:n_songs) {
  id_sx = annotations$song_id[sx];
  features_fname_sx = file.path(features_dir, paste0(id_sx, ".csv") );

  features_sx = read.csv(features_fname_sx, header=TRUE, sep=";");

  time_labels_map_sx = time_labels_map_init;

  # Use only the columns corresponding to times that also appear in the features table.
  ixs_time_to_keep = is.element(time_labels_map_sx$time, features_sx$frameTime);

  if (!any(ixs_time_to_keep) ) {
    warn_msg = "Features file times do not match the annotations file times for song %d. Skipping ..."
    warning(sprintf(warn_msg, id_sx) );

    next;
  }

  time_labels_map_sx = time_labels_map_sx[ixs_time_to_keep,];

  n_time_pts = nrow(time_labels_map_sx);

  # Stretch the annotations table.
  annotations_sx = pivot_longer(
    annotations[sx,c("song_id", time_labels_map_sx$time_label)],
    all_of(time_labels_map_sx$time_label),
    names_to="time_label",
    values_to="y");

  # Add the time.
  annotations_sx = merge(annotations_sx, time_labels_map_sx, by="time_label");

  # Filter the features table to just the time points for which there are responses present
  # in the annotations table and just the features needed to make the matrix-variate data.
  ixs_features = is.element(features_sx$frameTime, time_labels_map_sx$time);
  features_sx = features_sx[ixs_features,c("frameTime", feat_grp_labels)];
  rownames(features_sx) = NULL;

  # Initialize the map from the time in seconds to the corresponding matrix-variate features.
  features_map_sx = data.frame(time=time_labels_map_sx$time);
  features_map_sx$features = lapply(1:n_time_pts, \(ix) matrix(NA, n_feat_grps, len_feat_grps) );

  for (tx in 1:n_time_pts) {
    time_tx = features_map_sx$time[tx];
    ix_time = which(features_sx$frameTime == time_tx);

    features_tx = matrix(NA, n_feat_grps, len_feat_grps);

    for (fx in 1:n_feat_grps) {
      features_tx[fx,] = as.numeric(features_sx[ix_time,feat_grps[[fx]]]);
    }

    features_map_sx$features[[tx]][,] = features_tx;
  }

  # Match the annotations with the matrix-variate features by time.
  combined_datasets[[sx]] = merge(annotations_sx, features_map_sx, by="time");
}

combined_datasets = bind_rows(combined_datasets);

combined_datasets$time_bin = floor(combined_datasets$time);

# Thin the data by taking the mean across multiple time points.
thinned_data = (
  combined_datasets %>%
  group_by(song_id, time_bin) %>%
  summarise(
    y=mean(y),
    features=list(
      matrix(
        rowMeans(matrix(unlist(features), n_feat_grps * len_feat_grps) ),
        n_feat_grps) ),
    .groups="drop") );

# Correctness check.
#data_check1 = (
#  combined_datasets %>%
#  group_by(song_id, time_bin) %>%
#  nest() );

#ixs_order1 = order(data_check1$song_id, data_check1$time_bin);
#data_check1 = data_check1[ixs_order1,];

#ixs_order2 = order(thinned_data$song_id, thinned_data$time_bin);
#data_check2 = thinned_data[ixs_order2,];

#data_check_results = rep(NA, nrow(data_check1) );

#for (ix in 1:nrow(data_check1) ) {
#  res1 = data_check2$y[ix] == mean(data_check1$data[[ix]]$y);

#  accum_mat = matrix(0, n_feat_grps, len_feat_grps);

#  for (jx in 1:nrow(data_check1$data[[ix]]) ) {
#    accum_mat = accum_mat + data_check1$data[[ix]]$features[[jx]];
#  }

#  accum_mat = accum_mat / nrow(data_check1$data[[ix]]);
#  res2 = all.equal(data_check2$features[[ix]], accum_mat);

#  data_check_results[ix] = res1 && res2;
#}

#stopifnot(all(data_check_results) );

# Create covariates array.
X = array(
  unlist(thinned_data$features),
  dim=c(n_feat_grps, len_feat_grps, nrow(thinned_data) ));

# Correctness check.
#data_check_results = rep(NA, nrow(thinned_data) );

#for (ix in 1:nrow(thinned_data) ) {
#  data_check_results[ix] = all(thinned_data$features[[ix]] == X[,,ix]);
#}

#stopifnot(all(data_check_results) );

X_scaled = array(
  t(scale(t(matrix(X, n_feat_grps * len_feat_grps) ))),
  dim=dim(X) );

# The group_info table is an nx3 table where n is the total number of observations.
# The columns of the group_info table give the song ID, the time bin, and the response.
# The X array is p1xp2xn where p1 is the number of feature groups, an p2 is the length
# of the feature groups.
results = list();
results$group_info = select(thinned_data, song_id, time_bin, y);
results$X = X_scaled;

saveRDS(results, save_filename);
