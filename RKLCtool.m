function [adW, MW, W, AC, MC, PAC] = RKLCtool(X, T, kernel, ksigma)

% Do Regularized Kernel Logistic Classification 
% Input matrix X =>  Observation x  Features
% Input vector T =>  Label for each observation (Numeric values: 1~K)
%
%  AC: accuracy (%)
%  ac: pairwise accuracy
%
%        H.Oya ( modified , 2012) 
%% Do pairwise modeling 

K=length(unique(T));  % Number of class.........
iii=0;
W=zeros(size(X,2)+0,K*(K-1)/2)   ; R=[]; ii=0; Ti=[];
% Do all pairwise classifications .......
for n=1:K-1
    for m=n+1:K    
          Xi=[]; Ti=[];
          iii=iii+1;
          p=zeros(length(T),1);
         % ---------- Select stim. ID ------------ %
          r1=n; r2=m;
          fprintf('Classify %1.0f-%1.0f  \n ',r1,r2);
         %  ------------------------------------%
            f1=find(T==r1);  f2=find(T==r2); 
            Ti=T([f1; f2]);   
            f3=find(Ti==n);   f4=find(Ti==m);
            Ti(f3)=1; Ti(f4)=2;  % target vector is 0 or 1........
            Xi=X([f1; f2],:);
            XG=centernormalize(Xi,1);

            %  Find rambda -------------------------
            % using LOO accuracy (fast computing algorithm)
           [lambda] = findparameter_RKLR(XG, Ti, kernel, ksigma,[-2.5 0.5]);
           
            %  Model fitting ------------------------
           [alpha, beta, pp, w, Z, lik, ac(iii), Km ] = ...
                     Regu_Kernel_Logistic_Regress(XG, Ti, kernel, ksigma, lambda);
                 % beta is dual space weights. primal space
                 % weights is w
            %  Evaluate all points--------------------

%             W(:,iii)=[beta; w'];
            W(:,iii)=w';
%             eta=[ones(size(X,1),1) X]*W(:,iii); 
             eta=X*W(:,iii);
% size(X), size(alpha), size(beta)
%             eta = repmat(beta(1,:),size(X,1),1) + X*alpha;
            p=1./(1+exp(-eta));     % Posterior probabilty estimates..............
            p([f1; f2])=pp;         % Replace p by cross-validated one......
            R(n,m,:)=p;        % Make probability matrix.....
            clear pp p alpha beta w ;
    end
end

if K>2
    MW=mean(abs(W),2);
    adW=mean(abs(W).*repmat(ac(:)',size(W,1),1),2);
elseif K==2
    plot(squeeze(R(1,2,:)));
    MW=mean((W),2);
    f3=find(T==1);   f4=find(T==2);    
    ttt=length(f4)/(length(T));
    f1=find(R(1,2,:)<=0.5);    f2=find(R(1,2,:)>0.5);

    ac=(length(intersect(f1,f3))+length(intersect(f2,f4)))/length(T); 
    AC.ec=ac; AC.pc=ac; MC=[];PAC=[];
   adW=mean((W).*repmat(ac(:)',size(W,1),1),2);
end    


%% Pairwise-Coupling for multi-class classification
if K>2
A=[]; B=[]; r=[]; P=[]; 
for n=1:size(R,3)
     A=R(:,:,n);
     A=[A; zeros(1, size(A,2))];
     r=A+triu(1-A,1)';
     B=-r.*r';
     for k=1:size(B,2);
            B(k,k) =  sum(r(:,k).*r(:,k));
     end
     dimb=size(B,1);
     Z2=[B ones(dimb,1); ones(dimb,1)' 0];
     Z3=[zeros(dimb,1) ;1];
     P=Z2\Z3;
     P=P(1:dimb);
     [I,ii(n)]=max(P);
end

PC=100*length(find((ii'-T)==0))/length(T);
clear i A B Z2 Z3 dimb  m r
close all;

%%  Error-correction output coding for multi-class classificatoion
clear tem EC ;
EC=zeros(size(R,3),K*(K-1)/2);
fff=find(R==0.5);
lfff=length(fff);
R(fff)=R(fff)+randn(lfff,1)*10^-10;

for i=1:size(R,3)
e=0;
tem=[R(:,:,i) ;zeros(1, size(R,2))];
id=T(i);
    for n=1:K-1
        for m=n+1:K  
            e=e+1;
            EC(i,e)=-sign(tem(n,m)-0.5);   % Hamming
        end
    end
end

% Basecode-----------------
Basecode=zeros(K, K*(K-1)/2);
bf=0;
for n=1:K-1
    for m=n+1:K  
        bf=bf+1;
          Basecode(n,bf)=1;
          Basecode(m,bf)=-1;
    end
end

%  Hamming Distance decoding -------------------
Lb=size(Basecode,1);    D=[];
for n=1:size(EC,1)
     temp=repmat(EC(n,:),Lb,1).*Basecode;
     D(:,n)=sum((1-sign(temp))/2,2);
end

[VS, IS]=min(D,[],1);
ECOC=100*length(find((IS'-T)==0))/length(T);
clear n m e bf temp 

%% Confusion matrix 

du=[];si=[]; f=[]; 
for n=1:K
    f=find(T==n);L=length(f);
    du=length(si)+1: length(si)+L;
    si(du)=n;
end

 CL=length(find((si-IS)==0));
 AC.ec=CL*100/length(T);
%  CL=length(find((si-ii)==0));
 AC.pc=CL*100/length(T);
 
 fprintf(' Accuracy = %2.2f  percent \n', CL*100/length(T));
Oa=0; MC=[];MAC=[];
for n=1:K
    f1=find(si==n); g1=find(si~=n);
    for m=1:K
        % based on ECOC
        f2=find(IS(f1)==m);  g2=find(IS(g1)~=n);
        lf2=length(f2);
        MC(n,m)=lf2/length(f1); % true positive rate
        MAC(n,m)=(length(f2)+length(g2))/(length(T)); % Accuracy
        Oa=Oa+lf2;
    end
end

figure(24); imagesc(MC); set(gca,'ydir','norm','fontsize',14,'fontweight','bold')
load cm9; colormap(cm9);
if K~=2 
    caxis([1/K 0.45]);colorbar;
elseif k==2
    caxis([1/K 0.9]);colorbar;
end
ylabel(' Presented stim'); xlabel(' Predicted stim');
set(gcf,'paperunits','inches','papersize',[6 5],...
 'PaperPositionMode','manual','paperunits','normalized','paperposition',[0 0 1 1]);
clear tem Oa k f1 temp m f2 lf2 eta 


% eval(['print -f24 -dpng -r150 /home/hiroyuki/Tone/mitch/CM153ECOC-win',num2str(kk),';'])       

% figure(25); imagesc(MAC); set(gca,'ydir','norm','fontsize',14,'fontweight','bold')
% colormap(cm9);caxis([0 1]);colorbar;
% ylabel(' Presented stim'); xlabel(' Predicted stim');
% set(gcf,'paperunits','inches','papersize',[6 5],...
%  'PaperPositionMode','manual','paperunits','normalized','paperposition',[0 0 1 1]);
% eval(['print -f25 -dpng -r150 /home/hiroyuki/Tone/mitch/CM153PC-win',num2str(kk),';'])
% clear tem Oa k f1 temp m f2 lf2 eta  base 


%  Pairwise confusion matrix  /////////////////////////////
PAC=[];
i=0;
for n=1:K
    for m=n+1:K
        i=i+1;
        PAC(n,m)=ac(i);
    end
end
PAC=[PAC;zeros(1,K)];
tem=triu(PAC,1);
PAC=PAC+tem';

pPAC=triu(PAC);
figure(26); imagesc(pPAC/100); set(gca,'ydir','norm','fontsize',14,'fontweight','bold');
load cm10; 
colormap(cm10);caxis([0.55 0.95]);colorbar;title('Pairwise classification accuracy');
i=0; 
for n=1:K
    for m=n+1:K
        i=i+1;
        st=sprintf('%1.2f',ac(i));
        text(m-0.2,n-0.1 ,st,'color','b','fontsize',12,'fontweight','bold');
    end
end
set(gcf,'paperunits','inches','papersize',[6 5],...
 'PaperPositionMode','manual','paperunits','normalized','paperposition',[0 0 1 1]);
% eval(['print -f26 -dpng -r150 /home/hiroyuki/Tone/pac'])
end