function [W, P, AC] = HSLR_fit(X,T,L1, L2,Ns)
% functions needed:
%
%   1, Sparse_LR_CV.m
%   2, Sparse_MNL_Regression.m
%   L1=[a b]    L1 regularization parameter search region
%   L2=[c d]   L2 reguralization parameter search region
%   
%     H.Oya (2011)
%
%% =======  Parameter optimization =======%%
loomode=0;
ac=[];
if isempty(L1)==1
    o=-6;
else
    o=linspace(L1(1),   L1(2), Ns);     % L1 norm regularization
end

if  isempty(L2)==1
    o2=-6;
else
    o2=linspace(L2(1),  L2(2) ,Ns);    % L2 norm reguralization    
end
%==================================%
i=0;

Oind1=[];Oind2=[];tem=[];
L1=length(o);
L2=length(o2);
Oind2=repmat(o2,1,L2);
tem=o;  tem=repmat(tem,L2,1); 
Oind1=reshape(tem,L1*L2,1)';

%Grid search.................
ac=zeros(L2,L1); lik=zeros(L2,L1);
parfor m=1:L1*L2
     [W, ac(m), P, lik(m)] = Sparse_LR_CV(X,T,[10^Oind1(m)  10^Oind2(m)]);        
end
ac=reshape(ac, L1,L2);


% visualize accuracy ------
fprintf('\n');
figure(101); subplot(2,1,1); imagesc(ac); ylabel('L2'); xlabel('L1');colorbar;
                   ca=caxis;  mcs=max([45 ca(1)]); mcs2=max([72 ca(2)]); d=ca(2)-ca(1);
                   caxis([mcs mcs2]);
H=fspecial('gaussian',2,5);   ac2= imfilter(ac,H,'replicate');
figure(101); subplot(2,1,2);imagesc(ac2); ylabel('L2'); xlabel('L1');
                   caxis([mcs mcs2]); colorbar                  
drawnow;


[mm]=max(max(ac));
[u,uu]=find(ac==mm); %AC=ac(u(end),uu(end));
% Selected parameters........
gamma2=10^o2(u(end));
gamma1=10^o(uu(end));
% MH=mean(reshape(ac(1:8,1:8),8*8,1));
% MHM=max(reshape(ac(1:8,1:8),8*8,1));

% if MH<55 && MHM<65
%     gamma2=10^o2(1);
%     gamma1=10^o(1);
%     AC=ac(1,1);
% end

% Final W and P ........
W=[];
[W, AC, P] = Sparse_LR_CV(X,T,[gamma1 gamma2]);        
fprintf('     L1= %2.3f  :   L2= %2.3f  : ',gamma1,gamma2);

clear u uu m n i temac 
%%  LOO-CV ===========================
clear  p n 
if loomode==1;
    P=[];Wtemp=[];
    i=unique(T);
    K = length(i);     
      for n=1:length(T) 
          train=setdiff([1:length(T)],n);
          [prob,Wtemp(:,n), loglik]=Sparse_MNL_Regression(X(train,:), T(train), [gamma1 gamma2]);
          tprob=exp([1 X(n,:)]*Wtemp(:,n));
          if K~=2
              P(:,n)=tprob./repmat(sum(tprob,2),[1 2]);
          elseif K==2
              P(:,n)=tprob./(1+tprob);
          end
      end

      f=find(T==1);   f2=find(T==2);
      E=zeros(length(T),1); E(f)=1; E(f2)=-1;
      EE=sign(P(1,:)'-0.5);
      AC=length(find(E.*EE==1))/length(T);
      % Final model=================
       [tem1,W]=...
         Sparse_MNL_Regression(X, T, [gamma1 gamma2]);
end
  
%%
  f=find(W(2:end)~=0);
%   fprintf('     Selected %2.0f/%2.0f  : LOO-Accuracy  %2.2f  : N-fold accuracy  %2.2f \n', length(f), size(X,2),100*AC,mm);
  fprintf('Selected %2.0f/%2.0f  : N-fold accuracy  %2.2f \n', length(f), size(X,2), AC);
  fprintf('________________________ \n');
  clear f f2 E EE  n train tprob u uu prob 
