dreamer_filename = './DREAMER.mat';
save_filename = './processed_DREAMER_by_freq_bands.mat';

n_sec_per_video = 60;
max_freq = 50;

dreamer_data = load(dreamer_filename);
dreamer_data = dreamer_data.DREAMER;

fs = dreamer_data.EEG_SamplingRate;

max_freq = min(max_freq, fs / 2);

n_dim = dreamer_data.noOfSubjects;
m_dim = dreamer_data.noOfVideoSequences;

n_scores = 4;

channames = dreamer_data.EEG_Electrodes;
n_chans = length(channames);

n_fft = n_sec_per_video * fs;
freq_orig = reshape(0:(1 / n_sec_per_video):(fs / 2), [], 1);

freq_band_info = table( ...
  reshape(1:5, [], 1), ...
  {'Delta'; 'Theta'; 'Alpha'; 'Beta'; 'Gamma'}, ...
  [0.5; 4; 8; 12; 30], ...
  [4; 8; 12; 30; max_freq], ...
  'VariableNames', {'freq_group_id', 'freq_group', 'lb_freq', 'ub_freq'});

freq_group_ids = nan(length(freq_orig), 1);

for bx = 1:size(freq_band_info, 1)
  ixs_bx = ( ...
    (freq_band_info.lb_freq(bx) <= freq_orig) & ...
    (freq_orig < freq_band_info.ub_freq(bx) ));

  freq_group_ids(ixs_bx) = freq_band_info.freq_group_id(bx);
end

ixs_valid_freq = ~isnan(freq_group_ids);

freq_group_ids = freq_group_ids(ixs_valid_freq);
[freq_group_map, freq_groups] = findgroups(freq_group_ids);

n_freq = length(freq_groups);

eeg_stim_logpower = nan([n_freq, n_chans, m_dim, n_dim]);
scores = nan([m_dim, n_scores, n_dim]);

for ix = 1:n_dim
  eeg_ix = dreamer_data.Data{ix}.EEG;

  for mx = 1:m_dim
    % We consider only the last 'n_sec_per_video' seconds of each EEG recording.
    eeg_ix_mx = eeg_ix.stimuli{mx}((end - n_fft + 1):end,:);

    % Estimate PSD using Welch's method by dividing the time series into eight segments
    % with 50% percent overlap between segments. These are the default settings.
    [power_ix_mx, freq_check] = cpsd(eeg_ix_mx, eeg_ix_mx, [], [], n_fft, fs);

    assert(mean((freq_orig - freq_check).^2) < 1e-16);

    % Compute log of power and then average across frequency bins.
    logpower_ix_mx = log(power_ix_mx(ixs_valid_freq,:) );
    logpower_freq_groups = splitapply(@mean, logpower_ix_mx, freq_group_map);

    eeg_stim_logpower(:,:,mx,ix) = logpower_freq_groups;
  end

  scores(:,1,ix) = dreamer_data.Data{ix}.ScoreValence;
  scores(:,2,ix) = dreamer_data.Data{ix}.ScoreArousal;
  scores(:,3,ix) = dreamer_data.Data{ix}.ScoreDominance;
  scores(:,4,ix) = arrayfun(@(ix_row) sqrt(sum(scores(ix_row,1:3,ix).^2) ), 1:m_dim);
end

save_struct = struct();
save_struct.eeg_stim_logpower = eeg_stim_logpower;
save_struct.scores = scores;
save_struct.freq = reshape(freq_band_info.freq_group, 1, []);
save_struct.chans = channames;
save_struct.score_labels = {'Valence', 'Arousal', 'Dominance', 'Length'};
save_struct.age = cellfun(@(data_ix) str2double(data_ix.Age), dreamer_data.Data);
save_struct.gender = cellfun(@(data_ix) data_ix.Gender, dreamer_data.Data, 'UniformOutput', false);

save(save_filename, '-struct', 'save_struct');
