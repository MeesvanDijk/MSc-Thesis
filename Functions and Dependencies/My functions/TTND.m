function data_TT = TTND(X,l,eps);

%% Input
% X: a N*1 cell array for storing third-order tensor data
% l: rank truncation 
% output: 1. data_TT: TT decomposition of data,
%         2. Rank   : Rank corresponds to TT-decomposition.
%% Task : Getting compact representation of input data in TT form
N = length(X);

%% Decompose the input data with TT decomposition
data_TT=cell(N,1); % Save TT decomposition results
Rank = cell(N,1);
fprintf('Decomposing the input data with TT decomposition, please wait!\n');
for i=1:N
    %X_cwt = cwt_EEG(X{i}, 256);
    X_cwt = X{i};
    tt = tt_tensor(X_cwt,eps); % tt- decomposition of each cell tensor
    tt_round = round(tt,eps,l); % rounding of rank to l
    G = core2cell(tt_round); 
    data_TT{i}=cell(3,1);
    data_TT{i}{1}=G{1};
    data_TT{i}{2}=G{2};
    data_TT{i}{3}=G{3}; % G are the tt-cores
    clear X_cwt
end
end