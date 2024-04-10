clc
close
clear
%% Load data

patients = [1 3 8];                                                 
channels = linspace(1,18,18);


for patient = patients
    if patient == 1
        listing = struct2cell(dir('/Users/meesvandijk/Documents/MEES2/MEES/Systems & Control/Msc Thesis/EEG_MATLAB/Datasets/chb01')); 
        temp = listing(1, 4:45);
        filenames = strings(1,length(temp));
        for i = 1:length(temp)
            filenames(1,i) = string(temp{i});
        end
        
        chb01_all_data = cell(size(filenames));
        chb01_all_labels = cell(size(filenames));
        for i = 1:numel(filenames)
            filename =  filenames(i);
            filename_char =  char(filename);
            data = normalize(h5read(filename, '/'+ string(filename_char(1:8)) + '/eeg')); %Normalized
            chb01_all_data{i} = data(:,channels);
            chb01_all_labels{i} = h5read(filename, '/'+ string(filename_char(1:8)) + '/labels');
        end
    elseif patient == 2
        listing = struct2cell(dir('/Users/meesvandijk/Documents/MEES2/MEES/Systems & Control/Msc Thesis/EEG_MATLAB/Datasets/chb02'));
        temp = listing(1, 3:37);
        filenames = strings(1,length(temp));
        for i = 1:length(temp)
            filenames(1,i) = string(temp{i});
        end
        
        chb02_all_data = cell(size(filenames));
        chb02_all_labels = cell(size(filenames));
        for i = 1:numel(filenames)
            filename =  filenames(i);
            filename_char =  char(filename);
            data = normalize(h5read(filename, '/'+ string(filename_char(1:8)) + '/eeg')); %Normalized
            chb02_all_data{i} = data(:,channels);
            chb02_all_labels{i} = h5read(filename, '/'+ string(filename_char(1:8)) + '/labels');
        end
    elseif patient == 3
        listing = struct2cell(dir('/Users/meesvandijk/Documents/MEES2/MEES/Systems & Control/Msc Thesis/EEG_MATLAB/Datasets/chb03'));
        temp = listing(1, 3:40);
        filenames = strings(1,length(temp));
        for i = 1:length(temp)
            filenames(1,i) = string(temp{i});
        end
        
        chb03_all_data = cell(size(filenames));
        chb03_all_labels = cell(size(filenames));
        for i = 1:numel(filenames)
            filename =  filenames(i);
            filename_char =  char(filename);
            data = normalize(h5read(filename, '/'+ string(filename_char(1:8)) + '/eeg')); %Normalized
            chb03_all_data{i} = data(:,channels);
            chb03_all_labels{i} = h5read(filename, '/'+ string(filename_char(1:8)) + '/labels');
        end
    elseif patient == 8
        listing = struct2cell(dir('/Users/meesvandijk/Documents/MEES2/MEES/Systems & Control/Msc Thesis/EEG_MATLAB/Datasets/chb08'));
        temp = listing(1, 3:end);
        filenames = strings(1,length(temp));
        for i = 1:length(temp)
            filenames(1,i) = string(temp{i});
        end
        
        chb08_all_data = cell(size(filenames));
        chb08_all_labels = cell(size(filenames));
        for i = 1:numel(filenames)
            filename =  filenames(i);
            filename_char =  char(filename);
            data = normalize(h5read(filename, '/'+ string(filename_char(1:8)) + '/eeg')); %Normalized
            chb08_all_data{i} = data(:,channels);
            chb08_all_labels{i} = h5read(filename, '/'+ string(filename_char(1:8)) + '/labels');
        end
    end
    

end


%% Experiment LOSO
%Experiments setup
repetitions = 5;                                                   %The number of times the experiment is carried out in order to create confidence intervals
R_range = 5;                                                        %CP ranks included in the experiment
sample_patients = [1 3 8];                                                  %Select the patients for which to sample
metrics_cell_loso = cell(repetitions,length(sample_patients),R_range);                     %Set scores for SHTM
nr_samples = 1000;                                                  %Dummy value since LOSO ignores this value
window = 4;                                                         %time window of a single sample in seconds
overlap = 3;                                                        %overlap of the time windows
%%
loso = 1;                                                           %If LOSO crossvalidation is desired set this to 1
n_factor = 10;                                                       %Ratio of positive negative samples for 50/50 set this to 1
CWT = 0;                                                            %If CWT = 1 the output will be 3D tensors, if set to 0 the output will be 2D
g = [5:7];
C = [2:6];
flagker = 1;
kfold = 5;

store_C = zeros(R_range, length(sample_patients));
store_g = zeros(R_range, length(sample_patients));

%Experiment
for i = 1:repetitions
    for patient = 1:length(sample_patients)
        if sample_patients(patient) == 1
            [X, X_t, y, y_t, break_idx] = Data_samplerV2(chb01_all_data, chb01_all_labels, nr_samples, window, overlap, loso, n_factor, CWT);
        elseif sample_patients(patient) == 2
            [X, X_t, y, y_t, break_idx] = Data_samplerV2(chb02_all_data, chb02_all_labels, nr_samples, window, overlap, loso, n_factor, CWT);
        elseif sample_patients(patient) == 3
            [X, X_t, y, y_t, break_idx] = Data_samplerV2(chb03_all_data, chb03_all_labels, nr_samples, window, overlap, loso, n_factor, CWT);
        elseif sample_patients(patient) == 8
            [X, X_t, y, y_t, break_idx] = Data_samplerV2(chb08_all_data, chb08_all_labels, nr_samples, window, overlap, loso, n_factor, CWT);
        end
        div = length(break_idx)-1;
        X_SVD = SVND(X,R_range);
        X_SVD_t =  SVND(X_t,R_range);

        for R = 1:R_range
            X_SVD_R = cellfun(@(x) {x{1}(:,1:R), x{2}(:,1:R)}', X_SVD, 'UniformOutput', false);
            X_SVD_tR = cellfun(@(x) {x{1}(:,1:R), x{2}(:,1:R)}', X_SVD_t, 'UniformOutput', false);
            Order = length(X_SVD{1});
            [scores, preds, all_y_test, all_probs, copt, gopt] = LOSOV2(X_SVD_R, X_SVD_tR, y, y_t, break_idx, div, R, flagker, g, C, n_factor, kfold, i, store_C(R,patient), store_g(R, patient));
            store_C(R, patient) = copt;
            store_g(R, patient) = gopt;
            metrics_cell_loso{i,patient,R} = {scores, preds, all_y_test, all_probs};
            disp(scores)

    
        end
    end
    disp(i)
end

%%
file_name_loso = "/Users/meesvandijk/Documents/MEES2/MEES/Systems & Control/Msc Thesis/EEG_MATLAB/Results/E2_DuSK_LOSO.mat"
%% 
save(file_name_loso, 'metrics_cell_loso')




%clear chb01_all_data chb01_all_labels chb02_all_data chb02_all_labels chb03_all_data chb03_all_labels chb08_all_data chb08_all_labels    %Free up RAM


%% Experiment LOPO
%Experiments setup
repetitions = 5;                                                   %The number of times the experiment is carried out in order to create confidence intervals
R_range = 10;                                                       %CP ranks included in the experiment
sample_patients = [1];                                          %Select the patients for which to sample
metrics_cell_lopo = cell(repetitions,1, 3);                     %Set scores for SHTM
R_loop = [6,8,10]

nr_samples = 2000;                                                  %Dummy value since LOSO ignores this value
window = 8;                                                         %time window of a single sample in seconds
overlap = 6;                                                        %overlap of the time windows
loso = 0;                                                           %If LOSO crossvalidation is desired set this to 1
n_factor = 10;                                                       %Ratio of positive negative samples for 50/50 set this to 1
CWT = 0;                                                            %If CWT = 1 the output will be 3D tensors, if set to 0 the output will be 2D
g = [6];
C = [8];
flagker = 1;
kfold = 5;
c_opt = 0;
g_opt = 0;
div = 3;
store_C = zeros(3,div);
store_g = zeros(3,div);

%Experiment
for i = 1:repetitions
    clear X1 X_t1 y1 y1_t1 X3 X_t3 y3 y_t3 X8 X_t8 y8 y_t8 X1_CP X3_CP X8_CP X_t1_CP X_t3_CP X_t8_CP
    [X1, X_t1, y1, y_t1, break_idx1] = Data_samplerV2(chb01_all_data, chb01_all_labels, nr_samples, window, overlap, loso, n_factor, CWT);
    [X3, X_t3, y3, y_t3, break_idx3] = Data_samplerV2(chb03_all_data, chb03_all_labels, nr_samples, window, overlap, loso, n_factor, CWT);
    [X8, X_t8, y8, y_t8, break_idx8] = Data_samplerV2(chb08_all_data, chb08_all_labels, nr_samples, window, overlap, loso, n_factor, CWT);
    
    X1_SVD = SVND(X1,R_range);
    X3_SVD = SVND(X3,R_range);
    X8_SVD = SVND(X8,R_range);

    X_t1_SVD = SVND(X_t1,R_range);
    X_t3_SVD = SVND(X_t3,R_range);
    X_t8_SVD = SVND(X_t8,R_range);
    count = 0;
    for R = 6:2:10
        count = count + 1;
        X1_SVD_R = cellfun(@(x) {x{1}(:,1:R), x{2}(:,1:R)}', X1_SVD, 'UniformOutput', false);
        X3_SVD_R = cellfun(@(x) {x{1}(:,1:R), x{2}(:,1:R)}', X3_SVD, 'UniformOutput', false);
        X8_SVD_R = cellfun(@(x) {x{1}(:,1:R), x{2}(:,1:R)}', X8_SVD, 'UniformOutput', false);

        X_t1_SVD_R = cellfun(@(x) {x{1}(:,1:R), x{2}(:,1:R)}', X_t1_SVD, 'UniformOutput', false);
        X_t3_SVD_R = cellfun(@(x) {x{1}(:,1:R), x{2}(:,1:R)}', X_t3_SVD, 'UniformOutput', false);
        X_t8_SVD_R = cellfun(@(x) {x{1}(:,1:R), x{2}(:,1:R)}', X_t8_SVD, 'UniformOutput', false);

        X_SVD = {X1_SVD_R, X3_SVD_R, X8_SVD_R};
        X_SVD_t = {X_t1_SVD_R, X_t3_SVD_R, X_t8_SVD_R};

        y = {y1, y3, y8};
        y_t = {y_t1, y_t3, y_t8};
        div = length(y_t);
        Order = length(X_SVD{1}{1});
        % clear X1 X_t1 y1 y1_t1 X3 X_t3 y3 y_t3 X8 X_t8 y8 y_t8 X1_CP X3_CP X8_CP X_t1_CP X_t3_CP X_t8_CP
        [scores, preds, all_y_test, all_probs, copt, gopt] = LOPO(X_SVD, X_SVD_t, y, y_t, div, R, flagker, g, C, kfold, i, store_C, store_g);
        store_C(R,:) = copt;
        store_g(R,:) = gopt;
        metrics_cell_lopo{i,1,count} = {scores, preds, all_y_test, all_probs};
        disp(scores)


    end
    disp(i)
end


file_name_lopo = "/Users/meesvandijk/Documents/MEES2/MEES/Systems & Control/Msc Thesis/EEG_MATLAB/Results/E2_DuSK_LOPO2.mat"

save(file_name_lopo, 'metrics_cell_lopo')



%clear chb01_all_data chb01_all_labels chb02_all_data chb02_all_labels chb03_all_data chb03_all_labels chb08_all_data chb08_all_labels    %Free up RAM




%% plot results LOSO
close all
%load data
E1_data = load(file_name_loso);
E1_data = E1_data.metrics_cell_loso();
% F1-score plot LOSO
ticks = 1:R_range;
F1_tot_mat = zeros(repetitions, R_range);
F1_tot = 0;
cf_tot = zeros(2,2);
colors = ['r', 'g', 'b'];

for patient = 1:length(sample_patients)
    y = zeros(repetitions, R_range);
    for i = 1:repetitions
        for R = 1:R_range
            div = length(cell2mat(E1_data{i, patient,R}{1}));
            for j = 1:div
                [F1, ~, ~, cf_matrix] = metrics(E1_data{i, patient,R}{3}{j}, E1_data{i, patient,R}{2}{j});
                %cf_tot = cf_tot + cf_matrix;
                F1_tot_mat(i,R) = F1;
            end
            %[F1, recall, precision] = metrics2(cf_tot);
            %disp(cf_tot)
            %cf_tot = zeros(2,2);
            %F1_tot_mat(i,R) = F1;
            
        end

    end

    y = F1_tot_mat
    x = 1:R_range;
    N = size(y,1);                                      % Number of ‘Experiments’ In Data Set
    yMean = mean(y);                                    % Mean Of All Experiments At Each Value Of ‘x’
    ySEM = std(y)/sqrt(N);                              % Compute ‘Standard Error Of The Mean’ Of All Experiments At Each Value Of ‘x’
    CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    yCI95 = bsxfun(@times, ySEM, CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of ‘x’
    figure(1)
    plot(x, yMean,'-o', 'color',colors(patient))                                  % Plot Mean Of All Experiments
    hold on
    %plot(x, yCI95+yMean,'-r')                           % Plot 95% Confidence Intervals Of All Experiments
    patch([x, fliplr(x)], [yCI95(1,:)+yMean fliplr(yCI95(2,:)+yMean)] , colors(patient), 'EdgeColor','none', 'FaceAlpha',0.25)
    hold on
    grid
    hold on
end

figure(1)
xticks(ticks)
ylim([0 1])
xlabel('CP-Rank', 'Interpreter','Latex','Fontsize',15)
ylabel('F1-score', 'Interpreter','latex', 'Fontsize', 15)
title('F1-score DuSK LOSO', 'Interpreter', 'latex', 'Fontsize', 15)
legend('chb01', '', 'chb03', '', 'chb08', 'Interpreter', 'latex', 'Fontsize',15)


% Accuracy plot LOSO
colors = ['r', 'g', 'b'];
for patient = 1:length(sample_patients)
    y = zeros(repetitions, R_range)
    for i = 1:repetitions
        for R = 1:R_range
            div = length(cell2mat(E1_data{i, patient,R}{1}));
            y(i,R) = sum(cell2mat(E1_data{i, patient,R}{1})) / div;
    
        end
    end
    disp(y)
    x = 1:R_range;
    N = size(y,1);                                      % Number of ‘Experiments’ In Data Set
    yMean = mean(y);                                    % Mean Of All Experiments At Each Value Of ‘x’
    ySEM = std(y)/sqrt(N);                              % Compute ‘Standard Error Of The Mean’ Of All Experiments At Each Value Of ‘x’
    CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    yCI95 = bsxfun(@times, ySEM, CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of ‘x’
    figure(2)
    plot(x, yMean,'-o', 'color',colors(patient))                                  % Plot Mean Of All Experiments
    hold on
    %plot(x, yCI95+yMean,'-r')                           % Plot 95% Confidence Intervals Of All Experiments
    patch([x, fliplr(x)], [yCI95(1,:)+yMean fliplr(yCI95(2,:)+yMean)] , colors(patient), 'EdgeColor','none', 'FaceAlpha',0.25)
    hold on
    grid
    hold on
end

figure(2)
xticks(ticks)
ylim([0 100])
xlabel('CP-Rank', 'Interpreter','Latex','Fontsize',15)
ylabel('Accuracy [\%]', 'Interpreter','latex', 'Fontsize', 15)
title('Accuracy DuSK LOSO', 'Interpreter', 'latex', 'Fontsize', 15)
legend('chb01', '', 'chb03', '', 'chb08', 'Interpreter', 'latex', 'Fontsize',15)

%% plot results LOPO
close all
R_range = 5
repetitions = 5;
%load data
E1_data = load(file_name_lopo);
E1_data = E1_data.metrics_cell_lopo();
% F1-score plot LOPO
ticks = 1:R_range;
F1_tot_mat = zeros(repetitions, R_range);
F1_tot = 0;
cf_tot = zeros(2,2);
colors = ['r', 'g', 'b'];
sample_patients = [1 3 8]
for patient = 1:length(sample_patients)
    y = zeros(repetitions, R_range);
    for R = 1:R_range
        for i = 1:repetitions
            %div = length(cell2mat(E1_data{i, 1,R}{1}))
            [~, ~, ~, cf_matrix] = metrics(E1_data{i, 1,R}{3}{patient}, E1_data{i, 1,R}{2}{patient});
            [F1, recall, precision] = metrics2(cf_matrix);
            F1_tot_mat(i,R) = F1;
            
        end

    end

    y = F1_tot_mat
    x = 1:R_range;
    N = size(y,1);                                      % Number of ‘Experiments’ In Data Set
    yMean = mean(y);                                    % Mean Of All Experiments At Each Value Of ‘x’
    ySEM = std(y)/sqrt(N);                              % Compute ‘Standard Error Of The Mean’ Of All Experiments At Each Value Of ‘x’
    CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    yCI95 = bsxfun(@times, ySEM, CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of ‘x’
    figure(3)
    plot(x, yMean,'-o', 'color',colors(patient))                                  % Plot Mean Of All Experiments
    hold on
    %plot(x, yCI95+yMean,'-r')                           % Plot 95% Confidence Intervals Of All Experiments
    patch([x, fliplr(x)], [yCI95(1,:)+yMean fliplr(yCI95(2,:)+yMean)] , colors(patient), 'EdgeColor','none', 'FaceAlpha',0.25)
    hold on
    grid
    hold on
end

figure(3)
xticks(ticks)
ylim([0 1])
xlabel('CP-Rank', 'Interpreter','Latex','Fontsize',15)
ylabel('F1-score', 'Interpreter','latex', 'Fontsize', 15)
title('F1-score DuSK LOPO', 'Interpreter', 'latex', 'Fontsize', 15)
legend('chb01-03', '', 'chb01-08', '', 'chb03-08', 'Interpreter', 'latex', 'Fontsize',15)

% Accuracy plot LOPO

colors = ['r', 'g', 'b'];
for patient = 1:length(sample_patients)
    y = zeros(repetitions, R_range);
    for R = 1:R_range
        for i = 1:repetitions
            y(i,R) = E1_data{i, 1,R}{1}{patient}
    
        end
    end

    x = 1:R_range;
    N = size(y,1);                                      % Number of ‘Experiments’ In Data Set
    yMean = mean(y);                                    % Mean Of All Experiments At Each Value Of ‘x’
    ySEM = std(y)/sqrt(N);                              % Compute ‘Standard Error Of The Mean’ Of All Experiments At Each Value Of ‘x’
    CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    yCI95 = bsxfun(@times, ySEM, CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of ‘x’
    figure(4)
    plot(x, yMean,'-o', 'color',colors(patient))                                  % Plot Mean Of All Experiments
    hold on
    %plot(x, yCI95+yMean,'-r')                           % Plot 95% Confidence Intervals Of All Experiments
    patch([x, fliplr(x)], [yCI95(1,:)+yMean fliplr(yCI95(2,:)+yMean)] , colors(patient), 'EdgeColor','none', 'FaceAlpha',0.25)
    hold on
    grid
    hold on
end

figure(4)
xticks(ticks)
ylim([0 100])
xlabel('CP-Rank', 'Interpreter','Latex','Fontsize',15)
ylabel('Accuracy [\%]', 'Interpreter','latex', 'Fontsize', 15)
title('Accuracy DuSK LOPO', 'Interpreter', 'latex', 'Fontsize', 15)
legend('chb01-03', '', 'chb01-08', '', 'chb03-08', 'Interpreter', 'latex', 'Fontsize',15)