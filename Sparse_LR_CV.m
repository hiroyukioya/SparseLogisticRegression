function [W, ac, PR, loglik ,trial] = Sparse_LR_CV(X,T,param)
%
%   ac:  N-fCV accuracy;
%   T={1,2}
%% N-fold CV----------------------------------------
% Numkber of N-f CV
iter=2;
% Number of fold .......
N=5;    
xl=size(X,1);  
L=floor(xl/N); 
prob=[]; W=[]; newT=[] ; trial=[];
i=unique(T);
K = length(i);             %  K = number of class
if K==2
    KK=1;
else
    KK=K;
end
j=0;

[u,uu]=sort(T);
T=T(uu);
X=X(uu,:);

g=unique(T);
f1=length(find(T==g(1)))/length(T);
f2=length(find(T==g(2)))/length(T);

 [X] = centernormalize(X,0);

for m=1:iter
    % Shuffle trials ...............
    RA=randperm(xl);
    % Shuffled data sets ..............
    XE=X(RA,:);
    TE=T(RA);
    f3=find(TE(1:2*L)==1); f4=find(TE(2*L+1:4*L)==1);f5=find(TE([1:L  4*L+1:5*L])==1);

    while length(f3)/(2*L)>f1+0.1 || length(f3)/(2*L)<f1-0.1 || ...
            length(f4)/(2*L)>f1+0.1 || length(f4)/(2*L)<f1-0.1 || ... 
           length(f5)/(2*L)>f1+0.1 || length(f5)/(2*L)<f1-0.1
        % Re-shuffle trials ..............
        RA=randperm(xl);
        % Re-shuffled data sets .............
        XE=X(RA,:);
        TE=T(RA);
        f3=find(TE(1:2*L)==1); f4=find(TE(2*L+1:4*L)==1);f5=find(TE([1:L  4*L+1:5*L])==1);
%         fprintf('.');
    end

        % /////////////////     Loop over splits    ///////////////////%
        for k=1:N
            j=j+1;
            % Devide traials into test- and training trials ..........
            if k~=N
                 Te=(k-1)*L+1:L*k; 
                 Train=setdiff([1:xl],Te);
            elseif k==N
                 Te=(k-1)*L+1:xl;
                 Train=setdiff([1:xl],Te);
            end
            
            if K~=2
                % Fit the model..............
                [temprob,W(:,:,j), loglik] = ...
                    Sparse_MNL_Regression(XE(Train,:),TE(Train), [param(1) param(2)]);
                % Evaluate Test data points...........
                tprob=exp([ones(length(Te),1) XE(Te,:)]*W(:,:,j)); 
                prob= [prob ; tprob./repmat(sum(tprob,2),[1 KK])];
            elseif K==2
                % Fit the model..............
                [temprob,W(:,j), loglik] = ...
                    Sparse_MNL_Regression(XE(Train,:),TE(Train), [param(1) param(2)]);
                % Evaluate Test data points...........
                W(1,:)=0;
                tprob=exp([ones(length(Te),1) XE(Te,:)]*W(:,j)); 
                prob=[prob; tprob./(1+tprob)];
            end   
            trial=[trial  RA(Te)];
            newT=[newT; TE(Te)];
        end
end
if K~=2
    W=median(W,3);
elseif K==2
    W=median(W,2);
end

for n=1:size(X,1);
    f=find(trial==n);
    PR(n)=mean(prob(f));
end

f=find(newT(:,1)==1);
f2=find(newT(:,1)==2);
IND=zeros(length(T),1);   IND(f)=1; IND(f2)=-1; 
p=sign(prob(:,1)-0.5);
ac=100*length(find(IND.*p==1))/length(newT);

ID=zeros(length(T),1);   ID(f)=1; ID(f2)=0; 
loglik=0;       
% for n=1:length(newT)
%     loglik=loglik+ID(n)*log(prob(n))+(1-ID(n))*log(1-prob(n));
% end
%     loglik=loglik-param(1)*sum(abs(W(2:end)))-param(2)*sum(W(2:end).^2);

%%