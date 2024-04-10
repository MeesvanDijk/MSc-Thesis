clc
close
clear
%% Load data
%Change path to your own dataset folder

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
repetitions = 10;                                                           %The number of times the experiment is carried out in order to create confidence intervals
R_range = 5;                                                                %CP ranks included in the experiment
sample_patients = [1 3 8];                                                  %Select the patients for which to sample
metrics_cell_loso = cell(repetitions,length(sample_patients),R_range);      %Create empty cell for all information relevant for the performance meterics

nr_samples = 1000;                                                          %Dummy value since LOSO ignores this value
window = 6;                                                                 %time window of a single sample in seconds
overlap = 5;                                                                %overlap of the time windows
loso = 1;                                                                   %If LOSO crossvalidation is desired set this to 1
n_factor = 10;                                                              %Ratio of positive negative samples for 50/50 set this to 1
CWT = 1;                                                                    %If CWT = 1 the output will be 3D tensors, if set to 0 the output will be 2D
g = [0:8];                                                                  %Range of values for the kernel parameter "sigma"                                                          
C = [-8:8];                                                                 %Range of values for the regularization parameter "C"
flagker = 0;                                                                %If flagker=0 the SHTM (linear) kernel is used, when flagker=1 the DuSK-RBF kernel is used
kfold = 5;                                                                  %Amount of folds for determining the optimal "sigma" and "C" parameter on the training set

store_C = zeros(R_range, length(sample_patients));                          %Empty array for storing the optimal "C" parameter for every patient-rank combination
store_g = zeros(R_range, length(sample_patients));                          %Empty array for storing the optimal "sigma" parameter for every patient-rank combination

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
        for R = 1:R_range

            X_TT = TTND(X,R,epsilon);
            X_tTT = TTND(X_t,R, epsilon);
            dimension = [length(X_TT{1}{1}), length(X_TT{1}{2}), length(X_TT{1}{3})]
            X_CP = ttcptensor(X_TT,R,R,dimension,truncation);
            X_CP_t = ttcptensor(X_tTT, R, R, dimension, truncation);

            Rtot=R*R
            Order = length(X_CP{1});
            [scores, preds, all_y_test, all_probs, copt, gopt] = LOSOV2(X_CP, X_CP_t, y, y_t, break_idx, div, Rtot, flagker, g, C, n_factor, kfold, i, store_C(R,patient), store_g(R, patient));
            store_C(R, patient) = copt
            store_g(R, patient) = gopt
            metrics_cell_loso{i,patient,R} = {scores, preds, all_y_test, all_probs}
            disp(metrics_cell_loso{i, patient,R}{1})
            cf_tot = 0;
            for j = 1:div
                [~, ~, ~, cf_matrix] = metrics(all_y_test{j}, preds{j} );
                cf_tot = cf_tot + cf_matrix;
            end

            [F1, recall, precision] = metrics2(cf_tot)
    
        end
    end
    disp(i)
end

%%
file_name_loso = "/Users/meesvandijk/Documents/MEES2/MEES/Systems & Control/Msc Thesis/EEG_MATLAB/Results/E3_TTMMK_LOSO3_R.mat"
%%
save(file_name_loso, 'metrics_cell_loso')




%clear chb01_all_data chb01_all_labels chb02_all_data chb02_all_labels chb03_all_data chb03_all_labels chb08_all_data chb08_all_labels    %Free up RAM


%% Experiment LOPO
%Experiments setup
repetitions = 10;                                                           %The number of times the experiment is carried out in order to create confidence intervals
R_range = 5;                                                                %CP ranks included in the experiment
sample_patients = [1 3 8];                                                  %Select the patients for which to sample (Only relevant for plotting in the LOPO case)
metrics_cell_lopo = cell(repetitions,1, R_range);                           %Create empty cell for all information relevant for the performance metrics

nr_samples = 2000;                                                          %If value if lower than maximum amount a subset of positive instances is used.
window = 6;                                                                 %time window of a single sample in seconds
overlap = 5;                                                                %overlap of the time windows
loso = 0;                                                                   %If LOSO crossvalidation is desired set this to 1
n_factor = 10;                                                              %Ratio of positive negative samples for 50/50 set this to 1
CWT = 1;                                                                    %If CWT = 1 the output will be 3D tensors, if set to 0 the output will be 2D
g = [0:8];                                                                  %Range of values for the kernel parameter "sigma"                                                          
C = [-8:8];                                                                 %Range of values for the regularization parameter "C"
flagker = 1;                                                                %If flagker=0 the SHTM (linear) kernel is used, when flagker=1 the DuSK-RBF kernel is usedkfold = 5;
c_opt = 0;                                                                  %Initiate variable for storing an optimal c value
g_opt = 0;                                                                  %Initiate variable for storing an optimal sigma value
div = length(sample_patients);                                              %The amount of fold LOPO uses                                                             
store_C = zeros(R_range,div);                                               %Empty array for storing the optimal "C" parameter for every rank
store_g = zeros(R_range,div);                                               %Empty array for storing the optimal "sigma" parameter for every rank

%Experiment
for i = 1:repetitions
    clear X1 X_t1 y1 y1_t1 X3 X_t3 y3 y_t3 X8 X_t8 y8 y_t8 X1_CP X3_CP X8_CP X_t1_CP X_t3_CP X_t8_CP
    [X1, X_t1, y1, y_t1, break_idx1] = Data_samplerV2(chb01_all_data, chb01_all_labels, nr_samples, window, overlap, loso, n_factor, CWT);
    [X3, X_t3, y3, y_t3, break_idx3] = Data_samplerV2(chb03_all_data, chb03_all_labels, nr_samples, window, overlap, loso, n_factor, CWT);
    [X8, X_t8, y8, y_t8, break_idx8] = Data_samplerV2(chb08_all_data, chb08_all_labels, nr_samples, window, overlap, loso, n_factor, CWT);


    
    for R = 1:R_range
        X1_CP = CP3D(X1,R,1);
        X3_CP = CP3D(X3,R,1);
        X8_CP = CP3D(X8,R,1);

        X_t1_CP = CP3D(X_t1,R,1);
        X_t3_CP = CP3D(X_t3,R,1);
        X_t8_CP = CP3D(X_t8,R,1);

        X_CP = {X1_CP, X3_CP, X8_CP};
        X_CP_t = {X_t1_CP, X_t3_CP, X_t8_CP};

        y = {y1, y3, y8};
        y_t = {y_t1, y_t3, y_t8};
        div = length(y_t);
        Order = length(X_CP{1}{1})
        [scores, preds, all_y_test, all_probs, copt, gopt] = LOPO(X_CP, X_CP_t, y, y_t, div, R, flagker, g, C, kfold, i, store_C, store_g);
        store_C(R,:) = copt
        store_g(R,:) = gopt
        metrics_cell_lopo{i,1,R} = {scores, preds, all_y_test, all_probs}
        for j = 1:div
            [F1, ~, ~, cf_matrix] = metrics(all_y_test{j}, preds{j})
        end



    end
end

%%
file_name_lopo = "/Users/meesvandijk/Documents/MEES2/MEES/Systems & Control/Msc Thesis/EEG_MATLAB/Results/E1_SHTM_LOPO_R.mat"
%%
save(file_name_lopo, 'metrics_cell_lopo')



%clear chb01_all_data chb01_all_labels chb02_all_data chb02_all_labels chb03_all_data chb03_all_labels chb08_all_data chb08_all_labels    %Free up RAM




%% plot results LOSO
close all
R_range = 5;
repetitions =5;
sample_patients = [1 3 8];
%load data
E1_data = load(file_name_loso);
E1_data = E1_data.metrics_cell_loso();
% F1-score plot LOSO
ticks = 1:R_range;
F1_tot_mat = zeros(repetitions, R_range);
precision_tot_mat = zeros(repetitions, R_range);
recall_tot_mat = zeros(repetitions, R_range);
F1_tot = 0;
cf_tot = zeros(2,2);
colors = ['r', 'g', 'b'];

for patient = 1:length(sample_patients)
    y = zeros(repetitions, R_range);
    for i = 1:repetitions
        for R = 1:R_range
            div = length(cell2mat(E1_data{i, patient,R}{1}));
            for j = 1:div
                [~, ~, ~, cf_matrix] = metrics(E1_data{i, patient,R}{3}{j}, E1_data{i, patient,R}{2}{j});
                cf_tot = cf_tot + cf_matrix;
            end
            [F1, recall, precision] = metrics2(cf_tot);
            cf_tot = zeros(2,2);
            F1_tot_mat(i,R) = F1;
            precision_tot_mat(i,R) = precision;
            recall_tot_mat(i,R) = recall;
            
        end

    end
    
    recall_mean = mean(recall_tot_mat)
    precision_mean = mean(precision_tot_mat)
    y = F1_tot_mat;
    y = [y;y];
    x = 1:R_range;
    N = size(y,1);                                      % Number of ‘Experiments’ In Data Set
    yMean = mean(y)                                    % Mean Of All Experiments At Each Value Of ;x’
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
xlabel('TT-Rank', 'Interpreter','Latex','Fontsize',15)
ylabel('F1-score', 'Interpreter','latex', 'Fontsize', 15)
title('F1-score TT-MMK LOSO', 'Interpreter', 'latex', 'Fontsize', 15)
legend('chb01', '', 'chb03', '', 'chb08', 'Interpreter', 'latex', 'Fontsize',15)


% Accuracy plot LOSO
colors = ['r', 'g', 'b'];
for patient = 1:length(sample_patients)
    y = zeros(repetitions, R_range);
    for i = 1:repetitions
        for R = 1:R_range
            div = length(cell2mat(E1_data{i, patient,R}{1}));
            y(i,R) = sum(cell2mat(E1_data{i, patient,R}{1})) / div;
    
        end
    end
    x = 1:R_range;
    N = size(y,1);                                      % Number of ‘Experiments’ In Data Set
    yMean = mean(y)                                    % Mean Of All Experiments At Each Value Of ‘x’
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
xlabel('TT-Rank', 'Interpreter','Latex','Fontsize',15)
ylabel('Accuracy [\%]', 'Interpreter','latex', 'Fontsize', 15)
title('Accuracy TT-MMK LOSO', 'Interpreter', 'latex', 'Fontsize', 15)
legend('chb01', '', 'chb03', '', 'chb08', 'Interpreter', 'latex', 'Fontsize',15)

%% plot results LOPO
close all
%load data
repetitions = 2;
R_range = 5;
E1_data = load(file_name_lopo);
E1_data = E1_data.metrics_cell_lopo();
% F1-score plot LOPO
ticks = 1:R_range;
F1_tot_mat = zeros(repetitions, R_range);
precision_tot_mat = zeros(repetitions, R_range);
recall_tot_mat = zeros(repetitions, R_range);
F1_tot = 0;
cf_tot = zeros(2,2);
colors = ['r', 'g', 'b'];
sample_patients = [1 3 8];
for patient = 1:length(sample_patients)
    y = zeros(repetitions, R_range);
    for R = 1:R_range
        for i = 1:repetitions
            %div = length(cell2mat(E1_data{i, 1,R}{1}))
            [F1, recall, precision, cf_matrix] = metrics(E1_data{i, 1,R}{3}{patient}, E1_data{i, 1,R}{2}{patient});
            %cf_tot = cf_tot + cf_matrix
            if isnan(F1) == 1
                F1 = 0;
            end
            F1_tot_mat(i,R) = F1;
            precision_tot_mat(i,R) = precision;
            recall_tot_mat(i,R) = recall;
        end
        % [F1, recall, precision] = metrics2(cf_tot);
        % cf_tot = zeros(2,2);
        % F1_tot_mat(i,R) = F1;
    end
    
    recall_mean = mean(recall_tot_mat)
    precision_mean = mean(precision_tot_mat)
    y = F1_tot_mat;
    y = [y;y];
    x = 1:R_range;
    N = size(y,1);                                      % Number of ‘Experiments’ In Data Set
    yMean = mean(y)                                    % Mean Of All Experiments At Each Value Of ‘x’
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
xlabel('TT-Rank', 'Interpreter','Latex','Fontsize',15)
ylabel('F1-score', 'Interpreter','latex', 'Fontsize', 15)
title('F1-score TT-MMK LOPO', 'Interpreter', 'latex', 'Fontsize', 15)
legend('chb01-03', '', 'chb01-08', '', 'chb03-08', 'Interpreter', 'latex', 'Fontsize',15)

% Accuracy plot LOPO

colors = ['r', 'g', 'b'];
for patient = 1:length(sample_patients)
    y = zeros(repetitions, R_range);
    for R = 1:R_range
        for i = 1:repetitions
            y(i,R) = E1_data{i, 1,R}{1}{patient};
    
        end
    end
    y = [y;y];
    x = 1:R_range;
    N = size(y,1);                                      % Number of ‘Experiments’ In Data Set
    yMean = mean(y)                                     % Mean Of All Experiments At Each Value Of  x’
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
xlabel('TT-Rank', 'Interpreter','Latex','Fontsize',15)
ylabel('Accuracy [\%]', 'Interpreter','latex', 'Fontsize', 15)
title('Accuracy TT-MMK LOPO', 'Interpreter', 'latex', 'Fontsize', 15)
legend('chb01-03', '', 'chb01-08', '', 'chb03-08', 'Interpreter', 'latex', 'Fontsize',15)