function [scores, pred, all_training_kernels, all_test_kernels, all_probs] = Train_Test(Data_CP_train, Data_CP_test, y_train, y_test, R, Order, flagker, g1, C1)
%DUSK_MEESV2 Summary of this function goes here
%   Detailed explanation goes here

n=size(Data_CP_train,1);
all_training_kernels = cell(length(g1),1);
all_probs = cell(length(C1), length(g1));
Ktrain=zeros(length(n),length(n));
n_train =  n;
%tic
for k = 1:length(g1)
    %disp('Training kernel start')
    for i = 1:n
        for j = 1:i
            %Ktrain(i,j)= Ker_fTTCP(Data_CP_train{i,1},Data_CP_train{j,1},Order,2^g1(k), l, flagker);
            Ktrain(i,j)= Ker_DuSK2(Data_CP_train{i,1},Data_CP_train{j,1},Order,R,flagker,2^g1(k));
            if i~=j
                Ktrain(j,i)= Ktrain(i,j);
            end
        end
    end
    all_training_kernels{k,1} = [(1:n)',Ktrain];
    %disp('Training kernel finished')
end
%toc

n = size(Data_CP_test,1);
Ktest=zeros(n,n_train);
all_test_kernels = cell(length(g1),1);

for k = 1:length(g1)
    for i = 1:n 
        for j = 1:n_train
            %Ktest(i,j)= Ker_fTTCP(Data_CP_test{i,1},Data_CP_train{j,1},Order,2^g1(k), l, flagker);
            Ktest(i,j)= Ker_DuSK2(Data_CP_test{i,1},Data_CP_train{j,1},Order,R,flagker,2^g1(k));
        end
    end
    all_test_kernels{k,1} = [(1:n)', Ktest];

end


scores = zeros(length(C1),length(g1));

for i = 1:length(C1)
    for j = 1:length(g1)
        cmd=['-c ', num2str(2^C1(i)), ' -t ', num2str(4),' -q', '-b', num2str(1)];%,'-w1', num2str(2) ];
        model= svmtrain(y_train, all_training_kernels{j}, cmd);
        [pred,temp, prob] = svmpredict(y_test, all_test_kernels{j}, model,'-q');
        scores(i,j) = double(temp(1));
        all_probs{i,j} = prob;
        clear temp  model;
    end
end
%disp(scores)


end

