function [scores, preds, all_y_test, all_probs, c_opt, g_opt] = LOSOV2(X_CP, X_CP_t, y, y_t, break_idx, div, R, flagker, g1, C1, n_factor, kfold, i, store_C, store_g)
%LOSO Summary of this function goes here
%   Detailed explanation goes here



Order = length(X_CP{1});
scores = cell(1,div);
preds = cell(1,div);
all_probs = cell(1,div);
all_y_test = cell(1,div);
y_tests = zeros(1,div);
break_idx_neg = [0 break_idx(2:end)];
c_opt = 0;
g_opt = 0;

extra_break_idx = [0 break_idx(2:end)];
pos_neg = round(length(X_CP)/2);
X_pos = X_CP(pos_neg+1:end);
X_neg = X_CP(1:pos_neg);
y_pos = y(pos_neg+1:end);
y_neg = y(1:pos_neg);

if i == 1
    tic
    [C_optt, g_optt] = Kfold_Cg(X_CP, y, g1, C1, kfold, Order, flagker, R);
    toc
else
    C_optt = store_C;
    g_optt = store_g;
end

for k = 1:div

    ind_test_pos = (break_idx(k)+1):break_idx(k+1);
    ind_train_pos = setdiff( 1:length(X_pos), ind_test_pos);
    ind_test_neg = (break_idx_neg(k)+1):break_idx_neg(k+1);
    ind_train_neg = setdiff( 1:length(X_neg), ind_test_neg);
    y_test_extra = (extra_break_idx(k+1) - extra_break_idx(k)) * n_factor;
    X_test = [X_neg(ind_test_neg); X_pos(ind_test_pos); X_CP_t(1:y_test_extra)];
    test_indx = randperm(length(X_test));
    X_test =  X_test(test_indx);
    y_test = [y_neg(ind_test_neg); y_pos(ind_test_pos); y_t(1:y_test_extra) ];
    y_t =  y_t(y_test_extra+1:end);
    X_CP_t = X_CP_t(y_test_extra+1:end);
    y_test = y_test(test_indx);
    X_train = [X_neg(ind_train_neg); X_pos(ind_train_pos)];
    train_indx = randperm(length(X_train));
    X_train =  X_train(train_indx);
    y_train = [y_neg(ind_train_neg); y_pos(ind_train_pos)];
    y_train = y_train(train_indx);

    c_opt = C_optt;
    g_opt  = g_optt;
 
    [score, pred, ~, ~, probs] = Train_Test(X_train, X_test, y_train, y_test, R, Order, flagker, g_optt, C_optt);
    
    scores{1,k} = score;
    preds{1,k} = pred;
    all_probs{1,k} = probs;
    all_y_test{1,k} = y_test;
    
end

end