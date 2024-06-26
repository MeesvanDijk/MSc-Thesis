function ValueOfKernel = Ker_DuSK(X1,X2,Order,R,flagker,g);
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%%  DuSK kernel computing
values = zeros(R,R);
if nargin==3
 ValueOfKernel=X1'*X2;
else
    if flagker==0
        S=ones(R,R); 
        for k=1:Order
            S=S.*(X1{k,1}'*X2{k,1});           % X1 and X2 is Dk*R size
        end
        ValueOfKernel=sum(sum(S));
    else
        ValueOfKernel=0;
        gamma = 1/(2*g^2);
        for i=1:R
            for j=1:R
                exp_sum = 0;
                for k=1:Order
                    exp_sum = exp_sum -sum((X1{k,1}(:,i)-X2{k,1}(:,j)) .^ 2);
                end
                ValueOfKernel = ValueOfKernel + exp(exp_sum*gamma) ;
                
            end
            
        end
    end
end
end
