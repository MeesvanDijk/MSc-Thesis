function [X, X_test, y, y_test, break_idx_seizure] = Data_sampler(all_data, all_labels, nr_samples, window, overlap, loso, n_factor, CWT)
%DATA_SAMPLER Samples data points for a patient. The returned data can be
%be suitable for Leave One Seizure Out (LOSO) cross-validation is the loso
%input argument is set to 1.

%   Detailed explanation goes here
    rec_count = 0;
    fs = 256;
    window = window*fs; %convert seconds to total amount of samples
    overlap = overlap*fs;
    recording_amount = length(all_data);
    pos_class = cell(1,2000);
    neg_class = cell(1,20000);
    pos_idx = 0;
    neg_idx = 0;
    nr_neg_samp = 750;
    break_idx_seizure = zeros(1,40);
    if loso == 1
        nr_samples = 5000;
    end

    for i = 1:recording_amount
        recording = all_data{i};
        labels =  all_labels{i};
        dt = window - overlap;
        if any(labels > 0)
            seizure_start = find(labels > 0,1);
            seizure_duration = find(labels(seizure_start:end) < 0, 1) -1; %Works only if only one seizure in single recording
            
            
            steps = 1:floor((seizure_duration-dt)/dt);
            start_idx = seizure_start ;
            for j = 1:length(steps)
                end_idx = round(start_idx + window-1);
                sample = recording(start_idx:end_idx,:);
                tensor_sample = sample;
                pos_class{1,j+pos_idx} = tensor_sample;
                start_idx = start_idx + dt;
            end
            pos_idx = pos_idx + length(steps);
            break_idx_seizure(1,i) = pos_idx;
        else
            for j = 1:nr_neg_samp
                start_idx = max([round(rand * (length(recording))-window) 1]);
                end_idx = round(start_idx + window -1);
                sample = recording(start_idx:end_idx, :);
                tensor_sample = sample;
                neg_class{1,j+neg_idx} = tensor_sample;
            end
        neg_idx = neg_idx + nr_neg_samp;
        end
    rec_count = rec_count + 1;
    % clc
    % percent = round(rec_count/recording_amount*100);
    % mess1 = sprintf('Loading... %d%% ', percent);
    % disp(mess1);
    end
    
    break_idx_seizure = break_idx_seizure(find(break_idx_seizure));
    pos_class = pos_class(~cellfun('isempty',pos_class));
    neg_class = neg_class(~cellfun('isempty',neg_class));

    max_samples = length(pos_class) * 2;  %Define maximal amount of samples for training
    if max_samples >= nr_samples
        max_samples = round(nr_samples /2);
    else
        max_samples = max_samples / 2;
        mess2 = sprintf('Desired amount of samples is not possible. The maximal amount %d is used', max_samples);
        disp(mess2)
    end
    
    y1 = ones((max_samples),1)*-1;
    y2 = ones(max_samples,1);
    y_test = ones((n_factor*max_samples),1)*-1;
    y = [y1; y2];
    rand_indx = randperm(max_samples*(n_factor+1));
    neg_class_tot = neg_class(rand_indx);


    neg_class = neg_class_tot(1:max_samples);
    X_test = neg_class_tot(max_samples+1:end)';
    break_idx_seizure = [1 break_idx_seizure];

    if loso == 1 %If we need to do Leave One Seizure Out crossvalidation, we do not want the samples to be ordered randomly yet
        X = [neg_class'; pos_class'];
    else
        X = [neg_class'; pos_class'];
        break_idx = max_samples;
    end

end
