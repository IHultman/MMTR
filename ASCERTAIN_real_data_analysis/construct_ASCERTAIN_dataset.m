ascertain_features_fname = "./ASCERTAIN_Features/Dt_EEGFeatures.mat";
ascertain_ratings_fname = "./ASCERTAIN_Features/Dt_SelfReports.mat";

save_filename = "./ASCERTAIN_dataset.mat";

ascertain_features = load(ascertain_features_fname);
ascertain_ratings = load(ascertain_ratings_fname);

n_subj = 58;
n_video = 36;
n_features = 11;
n_chans = 8;

p_dim = n_features * n_chans;
n_tot_obs = n_video * n_subj;

assert(numel(ascertain_features.EEGFeatures_58) == n_subj);
assert( ...
  all( ...
    cellfun( ...
      @(feat_sx) size(feat_sx, 1) == n_video, ...
      ascertain_features.EEGFeatures_58) ));

assert( ...
  all( ...
    cellfun( ...
      @(feat_sx) size(feat_sx, 2) == p_dim, ...
      ascertain_features.EEGFeatures_58) ));

X = nan([p_dim, n_video, n_subj]);

for sx = 1:n_subj
  feat_sx = ascertain_features.EEGFeatures_58{sx}';
  X(:,:,sx) = feat_sx;
end

X = reshape(X, [p_dim, n_tot_obs]);
X(isinf(X) ) = nan;
X = normalize(X, 2);
X(isnan(X) ) = 0;
X = reshape(X, [n_features, n_chans, n_tot_obs]);

group_ids = repelem(1:n_subj, n_video);

y = reshape( ...
  permute(ascertain_ratings.Ratings, [1, 3, 2]), ...
  [5, n_tot_obs])';

y = normalize(y);
y(isnan(y) ) = 0;

save_struct = struct();
save_struct.X = X;
save_struct.group_ids = group_ids;
save_struct.y = y;

save(save_filename, '-struct', 'save_struct');
