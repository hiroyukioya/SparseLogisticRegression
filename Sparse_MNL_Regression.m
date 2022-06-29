function [prob,W, loglik]=Sparse_MNL_Regression(x,y,param)
%%  Sparse multinomial logistic regression 
%       < Hybrid (L1 and L2 ) regularization >
%
%  [ INPUTS ]:
%  X is an N-by-P design matrix:
%  Y is an N-by-1 vecot of response: (1~K)
%  param : regularization parametgers [lambda1 lambda2]
%                     lambda1=> L1 regularization value
%                     lambda2=> L2 regularization value
%   N: number of constraint (observation)
%   P: input feature dimension
%
%
%  [ Outputs ]:
% prob: predicted probability
% W:   Feature space weights 
%
%
%
%  [ References] : 
%
%         Krishnapurum  et.al.,    IEEE pattern anal & Machine intell (2005) 
%         Ryali et.al.,   Neuroimage (2010)
%         Le Cessie ;  J. R. Stat. Soc. C (1992)
%         G. Fort;       Bioinformatics (2005)
%         P. Green;   J. R. Stat. Soc. B (1984)
%         L. Shen;     IEEE/ACM. Comp. Biol. Bioinform. (2005)
%         J. Zhu;        Biostatistics (2004)
%         M. Park and Hastie;   Biostatistics (2008) 9,1, pp30-50
%         Cessie, S.;   Appl.Statist. (1992) 41, 191-201
%         Bull: 2002 ;   Comp stat and  Data analysis
%         Fox : Applied regression analysis and generalized linear models
%         (2008)
%         Bohning:   Ann.Inst. Statist. Math.  (1992) vol.44,  pp197-200
%         Cawlay, et.al,   Advances in neural information (2007) 

%          H.Oya (2010)

%%  Setting up the initial variables
prob=[];loglik=[];
X=[ones(size(x,1),1) x];    % Add constant term : X=[N x P];
y=y(:);
K=length(unique(y));  %  K = number of class
N=size(y,1);               %  N = number of observation
P = size(X,2);              %  P=dimention of features including intercept.
if K~=2
    Y = accumarray({(1:N)' y},1);    % Make Y matrix with indicators   
    W=0.1*ones(P,K);     % Initialize weights .....    w = [ P x  K]
    KK=K;
elseif K==2
    Y = accumarray({(1:N)' y},1);    % Make Y matrix with indicators   
    Y=Y(:,1);
    W=0.1*ones(P,1);     % Initialize weights .....    w = [ P x  K]
    KK=1;
end
oldWeights=W;
lambda1=param(1);   % regularization parameter (L1)
lambda2=param(2);   % regularization parameter (L2)

%%   Precompute B matrix  etc
B=X'*X;
%  We need only diagonal elements of this matrix.....
%  B = Lower bound of Hessian...
B=-0.5*((K-1)/K)*diag(B);      % B=[P x 1]  ( feature dim x 1) 

XY=X'*Y;          % [P x KK]     ( KK is number of Class)
Xw = X*W;       % [N x KK]
% E=1/K.*ones(N,KK);    % [N x KK]
E=exp(Xw);
if K~=2
    S=sum(E,2);    % [N x 1];
elseif K==2
        S=sum(E,2)+ones(N,1);    % [N x 1];
end
    
%  NOTE:   W = [Dimension x Class]  :

%% Main iteration {Component wise} ................................
iter=0;
tol=10^-3;
maxiter=1*10^4;

while iter<=maxiter
    iter=iter+1;
    for i=1:P   % Loop over feature dimension...
        for k=1:KK   % Loop over class... 
              oldW=oldWeights(i,k);
                if oldW~=0
                    pp=E(:,k)./S;
                    a=X(:,i)'*pp;
                    b=XY(i,k);

                    % Compute gradient...
                    grad=b-a;

                    % Update weights...
                    newW=(oldW*B(i)/(B(i)-lambda2))-(grad/(B(i)-lambda2));

                    % Compute thresholds ... 
                    thres1=-lambda1/(B(i)-lambda2);

                    % Aplly soft thresholding ....
                    newW=soft_threshold(newW, thres1);
                    
                    % Update E and S for next iteration.
                    if newW ~= oldW
                          Xw(:,k)=Xw(:,k)+X(:,i)*(newW-oldW);
                          Enew=exp(Xw(:,k));
                          S=S+(Enew-E(:,k));
                          E(:,k)=Enew;  
                          W(i,k)=newW;          
                    end
                    
%                     if newW-oldW)~=0
%                         W(i,k)=newW;  % Component wise update...
%                          tXw(:,k)=X(:,i)* newW;  % [N x1];
%                         E(:,k)=exp(Xw(:,k));      % [N x1]
%                         S=sum(E,2);
%                     end

            end
        end
    end
%     figure(100); plot(iter,log10(norm(W-oldWeights)/norm(oldWeights)),'.');hold on
   
    if norm(W-oldWeights)/norm(oldWeights)<tol
%            fprintf('%2.0f  iteration required.   \n',iter);
            if K~=2
                tprob=exp(X*W); 
                % posterior probabilities
                prob=tprob./repmat(sum(tprob,2),[1 KK]);
            elseif K==2
                prob=1./(1+exp(-X*W));
            end
            % likelihood
            loglik=sum(sum(log(prob).*Y));
        break
    elseif norm(W-oldWeights)/norm(oldWeights)>=tol  && iter<maxiter
        oldWeights=W;
    elseif norm(W-oldWeights)/norm(oldWeights)>=tol  && iter==maxiter
        disp(' iteration not converged......');
            if K~=2
                tprob=exp(X*W); 
                % posterior probabilities
                prob=tprob./repmat(sum(tprob,2),[1 KK]);
            elseif K==2
                prob=1./(1+exp(-X*W));
            end
            % likelihood
            loglik=sum(sum(log(prob).*Y));
        break
    end
end
        


        