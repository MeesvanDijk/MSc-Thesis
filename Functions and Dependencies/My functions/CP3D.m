function Data_CP = CP3D(X, R, init)
%If the CWT is used uncomment line 18 and comment line 17.


lsearch='elsr';     % line search, can be 'none', 'lsh', 'lsb', 'elsr' or 'elsc'
comp='on';          % ='on' or ='off' to perform dimensionality reduction or not
Tol1=1e-6;          % Tolerance of ALS 
MaxIt1=5e3;        % Max number of iterations 
Tol2=1e-4;          % tolerance in refinement stage (after decompression if it was used)
MaxIt2=5e1;         % Max number of iterations in refinement stage
Ninit=3;  

n=size(X,1);
Data_CP = cell(n,1);
if init==1
    for t=1:n
        X_cwt = cwt_EEG(X{t,1}, 256);
        %X_cwt = X{t,1};
        if t==1
            [A,B,C,~,~]=cp3_alsls(X_cwt,R,lsearch,comp,Tol1,MaxIt1,Tol2,MaxIt2,Ninit);
            Data_CP{t,1}=cell(3,1);
            Data_CP{t,1}{1,1}=A;
            Data_CP{t,1}{2,1}=B;
            Data_CP{t,1}{3,1}=C;
            A_init = A;
            B_init= B;
            C_init = C;
            clear X_cwt
        else
            [A,B,C,~,~]=cp3_alsls(X_cwt,R,lsearch,comp,Tol1,MaxIt1,Tol2,MaxIt2,Ninit, A_init, B_init, C_init);
            Data_CP{t,1}=cell(3,1);
            Data_CP{t,1}{1,1}=A;
            Data_CP{t,1}{2,1}=B;
            Data_CP{t,1}{3,1}=C;
            clear X_cwt
        end
    end

else
    for t=1:n
        %X_cwt = cwt_EEG(X{t,1}, 256);
        X_cwt = X{t,1};
        [A,B,C,~,~]=cp3_alsls(X_cwt,R);
        Data_CP{t,1}=cell(3,1);
        Data_CP{t,1}{1,1}=A;
        Data_CP{t,1}{2,1}=B;
        Data_CP{t,1}{3,1}=C;
        clear X_cwt
    end
end
clear A B C
end

