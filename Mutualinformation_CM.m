function [SSI, I, SP] = Mutualinformation_CM(MC,T)

% Mutual informaiton & SSI computation from confusion matrix 
% Compute:  I(S_pred; S)
% T:  stimulus code
% MC : Confusion matrix

for n=1:length(unique(T))
    N(n)=length(find(T==n));
end

pS=N./sum(N);
pRS=MC;

% Response Entropy
sm=sum(MC,1);
psm=sm./(sum(sm));
HR = -sum(psm.*log2(psm));

 % Stimulus  Entropy 
 HS = -sum(pS.*log2(pS));
 
% Noise entropy
for n=1:size(MC,1)
    f=find(pRS(n,:)~=0);
    nz(n)=length(f);
    HRS(n) = -sum(pRS(n,f).*log2(pRS(n,f)),2);
end
HRS=sum(pS.*HRS);

% Mutual information
I = HR - HRS;
%%
% Conditional p 
psr=MC./repmat(sm,size(MC,1),1);

for n=1:size(MC,1)
    f=find(psr(:,n)~=0);
    nzr(n)=length(f);
    HSR(n) = -sum(psr(f,n).*log2(psr(f,n)),1);
end

% Specific infiormation
SP=HS - HSR;

% Stimulus specific information
SSI = sum(pRS.*repmat(SP,size(MC,1),1),2);

