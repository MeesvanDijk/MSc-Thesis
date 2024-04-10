clc
close
clear
%% Load data

patients = [1];                                                 
channels = linspace(1,18,18);


for patient = patients
    if patient == 1
        listing = struct2cell(dir('EEG_MATLAB/Datasets/chb01')); 
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
        listing = struct2cell(dir('EEG_MATLAB/Datasets/chb02'));
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
        listing = struct2cell(dir('EEG_MATLAB/Datasets/chb03'));
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
        listing = struct2cell(dir('EEG_MATLAB/Datasets/chb08'));
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
R_range = 3;                                                        %CP ranks included in the experiment
sample_patients = [1 3 8];                                                  %Select the patients for which to sample
metrics_cell_loso = cell(repetitions,length(sample_patients),R_range);                     %Set scores for SHTM

nr_samples = 1000;                                                  %Dummy value since LOSO ignores this value
window = 4;                                                         %time window of a single sample in seconds
overlap = 1;                                                        %overlap of the time windows
loso = 1;                                                           %If LOSO crossvalidation is desired set this to 1
n_factor = 10;                                                       %Ratio of positive negative samples for 50/50 set this to 1
CWT = 1;                                                            %If CWT = 1 the output will be 3D tensors, if set to 0 the output will be 2D
g = [6];
C = [0];
flagker = 1;
kfold = 5;



%%
%R=10
%R_range = 10
[X, X_t, y, y_t, break_idx] = Data_samplerV2(chb01_all_data, chb01_all_labels, nr_samples, window, overlap, loso, n_factor, CWT);
%%
% X_SVD = SVND(X,R_range);
% X_SVD_t =  SVND(X_t,R_range);
% Order = length(X_SVD{1});
%div = length(break_idx)-1;
%%
n = 1000;
X = cell(n,1);
for i=1:n
    X{i,1} = rand(5,5,5,5,5);
end

R=10
X_CP = CPND(X,R);
%X_CP_t =  CP3D(X_t, R,1);
Order = length(X_CP{1});
%%
n=size(X_CP,1);
Ktrain_old=zeros(length(n),length(n));
Ktrain_new=zeros(length(n),length(n));
n_train =  n;

ranks = 10;
times_old = zeros(1,10)
times_new = zeros(1,10)
g = 6
C = 6
for R=1:10
    tic
    %X_CP = cellfun(@(x) {x{1}(:,1:R), x{2}(:,1:R)}', X_SVD, 'UniformOutput', false);
    %X_SVD_tR = cellfun(@(x) {x{1}(:,1:R), x{2}(:,1:R)}', X_SVD_t, 'UniformOutput', false);
    for i = 1:n
        for j = 1:i
            %Ktrain(i,j)= Ker_fTTCP(Data_CP_train{i,1},Data_CP_train{j,1},Order,2^g1(k), l, flagker);
            Ktrain_old(i,j)= Ker_DuSK(X_CP{i,1},X_CP{j,1},Order,R,flagker,2^g);
            if i~=j
                Ktrain_old(j,i)= Ktrain_old(i,j);
            end
        end
    end
    times_old(1,R) = toc
    
    tic
    for i = 1:n
        for j = 1:i
            %Ktrain(i,j)= Ker_fTTCP(Data_CP_train{i,1},Data_CP_train{j,1},Order,2^g1(k), l, flagker);
            Ktrain_new(i,j)= Ker_DuSK2(X_CP{i,1},X_CP{j,1},Order,R,flagker,2^g);
            if i~=j
                Ktrain_new(j,i)= Ktrain_new(i,j);
            end
        end
    end
    %Ktrain_new = Ktrain_new + tril(Ktrain_new,-1)';
    times_new(1,R) = toc
end
%% check if kernels are the same

check =  norm(Ktrain_old - Ktrain_new)
%K_dif = Ktrain_old - Ktrain_new
%%
close all
plot_ranks = 1:10
figure(1)
plot(plot_ranks, times_old)
hold on
plot(plot_ranks, times_new)
legend('Old method', 'New method', 'Interpreter', 'Latex', 'FontSize',15)
title('Kernel computation time for different CP-ranks','Interpreter', 'Latex', 'FontSize',15)
xlabel('CP-rank', 'Interpreter', 'Latex', 'FontSize',15)
ylabel('Time [s]', 'Interpreter', 'Latex', 'FontSize',15)

figure(2)
plot(plot_ranks, (times_old - times_new)./ times_old *100,'bo-')
title('Change in kernel computation time for different CP-ranks','Interpreter', 'Latex', 'FontSize',15)
xlabel('CP-rank', 'Interpreter', 'Latex', 'FontSize',15)
ylabel('$\Delta$ kernel computation time [$\%$]', 'Interpreter', 'Latex', 'FontSize',15)


%%

vec1 = rand(1,10000);
vec2 = rand(1,10000);

tic
norm(vec1-vec2)^2
toc

tic
sum((vec1-vec2).^2)
toc

%% SVD TEST
X = rand(5,5)

[U, S, V] = svd(X)
X_rec = U*S*V'

X_rec_sum = zeros(5,5)
for i=1:length(X)
    X_rec_sum = X_rec_sum + (U(:,i)*S(i,i)*V(:,i)')
end

check1 = norm(X-X_rec )
check2 = norm(X-X_rec_sum)