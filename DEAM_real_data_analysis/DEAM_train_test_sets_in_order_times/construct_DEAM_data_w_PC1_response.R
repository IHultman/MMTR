DEAM_arousal_filename = "./DEAM_preprocessed_data_averaged_across_rater_arousal_response.RData";
DEAM_valence_filename = "./DEAM_preprocessed_data_averaged_across_rater_valence_response.RData";

DEAM_save_filename = "./DEAM_preprocessed_data_averaged_across_rater_PC1_response.RData";

DEAM_arousal_data = readRDS(DEAM_arousal_filename);
DEAM_valence_data = readRDS(DEAM_valence_filename);

stopifnot(
  all(
    DEAM_arousal_data$group_info[,c("song_id", "time_bin")] ==
    DEAM_valence_data$group_info[,c("song_id", "time_bin")]) );

stopifnot(all(DEAM_arousal_data$X == DEAM_valence_data$X) );

response_mat = cbind(DEAM_arousal_data$group_info$y, DEAM_valence_data$group_info$y);
response_pca = prcomp(response_mat, center=TRUE, scale.=TRUE);

PC1_group_info = cbind(
  DEAM_arousal_data$group_info[,c("song_id", "time_bin")],
  y=response_pca$x[,1]);

DEAM_PC1_data = list(
  X=DEAM_arousal_data$X,
  group_info=PC1_group_info);

saveRDS(DEAM_PC1_data, DEAM_save_filename);
