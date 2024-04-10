function [C_opt, g_opt] = kfold_Cg(X_CP, y, g, C, kfold, Order, flagker, R)


n=size(X_CP,1);                                                                                            %Total number of samples
Kernel=zeros(length(n),length(n));                                                                         %Initiate empty kernel
all_kernels = cell(length(g),1);                                                                           %Save all kernels for all g/C combinations
all_acc = zeros(length(C), length(g));



for k = 1:length(g)
    for i = 1:n
        for j = 1:i
            Kernel(i,j)= Ker_DuSK(X_CP{i,1},X_CP{j,1},Order,R,flagker,2^g(k));
            if i~=j
                Kernel(j,i)= Kernel(i,j);
            end
        end
    end
    all_kernels{k,1} = [(1:n)',Kernel];
end

for i = 1:length(C)
    for j = 1:length(g)
        cmd=['-c ', num2str(2^C(i)), ' -t ', num2str(4), ' -v ' , num2str(kfold), ' -q '];                  % -v option is for the amount of folds
        Acc = svmtrain(y, all_kernels{j}, cmd);
        all_acc(i,j) = Acc;
    end
end

max_acc = max(all_acc, [], 'all');
[row, col] =  find(all_acc == max_acc,1);
C_opt = C(row);
g_opt = g(col);

end

