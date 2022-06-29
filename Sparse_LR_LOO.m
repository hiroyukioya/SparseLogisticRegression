close all;ac=[]
o=linspace(-2,.9,20);
lambda2 = 10^2;
h=waitbar(0, 'wait....');
for n=1:20 
    waitbar(n/20);
    [W, ac(n),p] = Sparse_LR_CV(X,T+1,[10^o(n) lambda2]);
end
close(h);
plot(10.^o,ac);set(gca,'xscale','log');grid on;
[u,uu]=max(ac);u

% final model
%  [W, ac,p] = Sparse_LR_CV(X,T+1,[10^o(11) 100]);
%   [prob,W, loglik]=Sparse_MNL_Regression(X, T+1, [ 10^o(uu) 100]);
  
%% =======  Parameter optimization ========
close all;ac=[];
o=linspace(-1.0,1.0,10);  % L1 norm regularization
o2=linspace(-2, 3,15);     % L2 norm reguralization

h=waitbar(0, 'wait....');
i=0;
for m=1:15
    for n=1:10
        i=i+1;
        waitbar(i/(15*10));
        [W, ac(m,n)] = Sparse_LR_CV(X,T+1,[10^o(n) 10^o2(m)]);
    end
end
close(h);
imagesc(ac); ylabel('L2'); xlabel('L1');colorbar


mm=max(max(ac)),
[u,uu]=find(ac==mm);
% [I,J]=ind2sub([15 15],uu);
gamma2=10^o2(u(end))
gamma1=10^o(uu(end))

  %% LOO-CV
clear W p n 
P=[];W=[];
i=unique(T);
K = length(i);     
  for n=1:length(T) 
      train=setdiff([1:length(T)],n);
      [prob,W, loglik]=Sparse_MNL_Regression(X(train,:), T(train)+1, [gamma1 gamma2]);
      tprob=exp([1 X(n,:)]*W);
      if K~=2
         P(:,n)=tprob./repmat(sum(tprob,2),[1 2]);
      elseif K==2
          P(:,n)=tprob./(1+tprob);
      end
  end
  f=find(T==0);
  f2=find(T==1);
  E=zeros(length(T),1); E(f)=1; E(f2)=-1;
  EE=sign(P(1,:)'-0.5);
  AC=length(find(E.*EE==1))/length(T);
  f=find(W(:,1)~=0);
  fprintf('Selected %2.0f out of %2.0f  : LOO-Accuracy  %2.2f  \n', length(f), size(X,2),AC);
  
  clear f f2 E EE  n train tprob u uu prob 