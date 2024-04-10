function [scores, preds, all_y_test, all_probs, c_opt, g_opt] = LOPO(X_CP, X_CP_t, y, y_t, div, R_tot, flagker, g1, C1, kfold, i, store_C, store_g)
%LOSO Summary of this function goes here
%   Detailed explanation goes here


R = round(sqrt(R_tot));
Order = length(X_CP{1}{1});
scores = cell(1,div);
preds = cell(1,div);
all_probs = cell(1,div);
y_tests = zeros(1,div);
c_opt = zeros(1,div);
g_opt = zeros(1,div);

patients =  1:length(X_CP);
combs = nchoosek(patients, (div-1))

for k = 1:div
    testp = setdiff(patients,combs(k,:));
    combs(k,1)
    combs(k,2)
    X_train = [X_CP{combs(k,1)}; X_CP{combs(k,2)}];
    X_test = [X_CP{testp}; X_CP_t{testp}];
    y_train = [y{combs(k,1)}; y{combs(k,2)}];
    y_test = [y{testp}; y_t{testp}];

    %Randomize order of samples and labels
    rndm_train = randperm(length(X_train));
    rndm_test = randperm(length(X_test));

    X_train = X_train(rndm_train);
    y_train = y_train(rndm_train);
    X_test = X_test(rndm_test);
    y_test = y_test(rndm_test);
    
    if i == 1
        tic
        [C_optt, g_optt] = Kfold_Cg(X_train, y_train, g1, C1, kfold, Order, flagker, R_tot);
        toc
        c_opt(1,k) = C_optt;
        g_opt(1,k) = g_optt;
    else
        C_optt = store_C(R,k);
        g_optt = store_g(R,k);
        c_opt(1,k) = C_optt;
        g_opt(1,k) = g_optt;
    end
    
    [score, pred, ~, ~, probs] = Train_Test(X_train, X_test, y_train, y_test, R_tot, Order, flagker, g_optt, C_optt);

    scores{1,k} = score;
    preds{1,k} = pred;
    all_probs{1,k} = probs;
    all_y_test{1,k} = y_test;
    
end

end

