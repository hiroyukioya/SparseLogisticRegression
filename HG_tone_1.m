%%  Do classification ...........
cd /home/hiroyuki/SLCtool;
X=[]; T=[];close all;
clear mw SSI I A;
eval(['mkdir ~/Tone/HG_tone/',pid,';']);
ch=[1:14];
intens=2;
for twin=1:length(TT)
%     twin=20;
    X=[]; T=[];
%     temX=squeeze(G(:,[1:4:24], twin, ch));
%     temX=squeeze(G(:, [1 5 17 25 41 45], twin, ch)); % 48 stim
     temX=squeeze(G(:, [intens:4:24], twin, ch));  % 24 stim
%      temX=squeeze(G(:, [3:4:24], twin, ch));  % 24 stim

    [a,b,c]=size(temX);
    X=reshape(temX,a*b,c);

    for n=1:size(temX,2) 
         T=[T;ones(size(temX,1),1)*n];
    end
    fi=find(isnan(sum(X,2))==1);
    T(fi)=[];
    X(fi,:)=[];

    [mw(:,twin), rw(:,twin), temp, ac, MC(:,:,twin), PMC(:,:,twin)] = SLCtool(X, T);
    [SSI(twin,:), I(twin)] = Mutualinformation_CM(MC(:,:,twin),T);
    A(twin)=ac.ec;
    set(gcf,'paperunits','inches','papersize',[6 5],...
    'PaperPositionMode','manual','paperunits','normalized','paperposition',[0 0 1 1]);
    eval(['print -f24 -dpng -r300 ~/Tone/HG_tone/',pid,'/',pid,'_cm_tw',num2str(twin),'_int',num2str(intens),';'])
    eval(['print -f26 -dpng -r300 ~/Tone/HG_tone/',pid,'/',pid,'_pcm_tw',num2str(twin),'_int',num2str(intens),';'])
    close all;
end
AW = repmat(A/100,size(mw,1),1).*mw; % accuracy adjusted weights

clear hdata s bi n  X twin rej* reft ss temX tempf tempw tstamp maxtri a Zin ... 
    ans c bias d fi f XM XI Ts Jr2 Jr Li YI ZI S data ch b cho g i j tem

eval(['save ~/Tone/HG_tone/',pid,'/',pid,'_classifier_int',num2str(intens),' ;']) 

%
cd ~/sfn2010
    [ZY,x,y] = Interp2d(SSI, [TT(1) TT(size(SSI,1))],  [1 6] , 150);
    figure(76);h1=subplot(1,2,1);imagesc(x,y,ZY');caxis([0.05 1.25]); 
    set(gca,'fontsize',15,'fontweight','bold')
    set(gcf,'position',[50 50 700 800])
    set(h1,'position',[0.14 0.1 0.4 0.82]); ylim([-0.1 0.5]);
    ylabel('sec'); xlabel('stim freq');
    
    h2=subplot(1,2,2); plot(I,TT,'linewidth',2);grid on;ylim([-0.1 0.5]);axis on
    set(h2,'position',[0.6 0.1 0.3 0.82],'yaxislocation','right','ydir','rev');
    xlabel('MI_lower bound (bits)');    set(gca,'fontsize',15,'fontweight','bold')
     set(gcf,'position',[50 50 720 850])
     
     set(gcf,'paperunits','inches','papersize',[6 7],...
 'PaperPositionMode','manual','paperunits','normalized','paperposition',[0 0 1 1]);

     eval(['print  -f76 -dpng -r150  ~/Tone/HG_tone/decoded_info_',pid,'_int',num2str(intens),' ;'])
        cd /home/hiroyuki/SLCtool;
