%%
clear
PID=[130;139;146;151;153;154; 138;162];
pat=8;
% for pat=7
    x=[]; T=[]; clear IW AW MW PAC MC block pid blk 
    pid=PID(pat);
    
    if pid==130
        block='061'; blk=1; reso=1;
    elseif pid==139
        block='019'; blk=1; reso=0;
    elseif pid==146
        block='004'; blk=1; reso=0;
    elseif pid==151
        block='035'; blk=1; reso=0;
    elseif pid==153
        block='019'; blk=1; reso=0;
    elseif pid==154
        block='005'; blk=1; reso=1;
    elseif pid==138
        block='004'; blk=1; reso=0;
    elseif pid==162
        block='019'; blk=1; reso=0;     
    end
      pid=num2str(pid);
    
    for w=1:16
          ww=w;
          eval(['nx',num2str(ww),'=[];']);
          T=[];
        for n=1:6
            i=n-1;
            if reso==1
            eval(['a1=importdata(''D:\Tone\GT8s\All_trials_ERBP_',pid,'-',block,'_blk',num2str(blk),'_',num2str(i),'_trials_resort_win',num2str(w),'.txt'');'])
            else
            eval(['a1=importdata(''D:\Tone\GT8s\All_trials_ERBP_',pid,'-',block,'_blk',num2str(blk),'_',num2str(i),'_trials_win',num2str(w),'.txt'');'])
            end
            sx1x=a1(:,[2:end])';    ch=a1(:,1);
            clear a1; [a,b]=size(sx1x);
            eval(['nx',num2str(ww),'=[nx',num2str(ww),';sx1x];']);   clear sx1x;
            tempT=ones(a,1)*n;  T=[T; tempT];
        end     
        eval(['x',num2str(ww),'=centernormalize(nx',num2str(ww),',1);']);
    end

    twin=[-0.1:0.05:0.65]+0.025;
    clear  a temp i n s s2 tempT b n w ww i temp* x 
    eval(['save D:\Tone\GridToneAnal_pt',pid,'  T x* ch nx* twin;']);
%%
    %%%%%%%%%%%%%%%%%%%%
    %%%%%%   LOAD DATA  %%%%%%%
    %%%%%%%%%%%%%%%%%%%%

    eval(['load  D:\Tone\GridToneAnal_pt',pid,' ;'])
    cd D:\SLCtool;
    eval(['mkdir D:\Tone\',pid,';'])
    
    %%%%%%  Prefilter features  %%%%%%%
    for n=1:16
         eval(['Mm(:,n)=mean(nx',num2str(n),',1);'])
    end
    base=mean(Mm(:,1:3),2);
    Mm=Mm-repmat(base,[ 1 size(Mm,2)]);
    figure(22); plot(Mm');grid on;
    selch=find(max(abs(Mm(:,1:10)),[],2)>0.5);
    sprintf('selected features = % 2.1d',length(selch))
    clear base Mm
    
    %% ===========================%
    %%%%%%  Do Classification:  %%%%%%%
    clear nx*; close all;    
%     matlabpool close;
    matlabpool open 4
    T=T(:);
    K=length(unique(T));     % Number of class
    MW=zeros(size(x1,2),16);      AW=zeros(size(x1,2),16);
    IW=zeros(size(x1,2),K*(K-1)/2, 16);   MC=zeros(K,K,16);
    PAC=zeros(K,K,16);
    aIW=zeros(size(x1,2),K*(K-1)/2, 16);  aW=zeros(size(x1,2),16);
    
    hh=waitbar(0, ' Calculating ....');
    for tw=1:16
        waitbar(tw/16);
        eval(['X=x',num2str(tw),'(:,selch);']); 
        [adw, mw, tempw, AC, MC(:,:,tw), PAC(:,:,tw)] = SLCtool(X, T);
        eval(['print -f24 -depsc -r150 D:\Tone\',pid,'\',pid,'_cm_tw',num2str(tw),';'])
        eval(['print -f26 -depsc -r150 D:\Tone\',pid,'\',pid,'_pcm_tw',num2str(tw),';'])
        eval(['print -f24 -dpng -r300 D:\Tone\',pid,'\',pid,'_cm_tw',num2str(tw),';'])
        eval(['print -f26 -dpng -r300 D:\Tone\',pid,'\',pid,'_pcm_tw',num2str(tw),';'])

        % Recover original W 
        MW(selch,tw)=mw(2:end,:);   
        AW(selch,tw)=adw(2:end,:);   
        IW(selch,:, tw)=tempw(2:end,:);
        Accuracy(:,tw)=[AC.ec; AC.pc];
        aIW(:,:,tw)=0.01*Accuracy(1,tw).*IW(:,:,tw);
        aW(:,tw)=0.01*Accuracy(1,tw).*MW(:,tw);
        close figure 24; close  figure 26
    end
    close all;
    matlabpool close
    eval(['save D:\Tone\',pid,'\',pid,'_classification_result MW AW  IW Accuracy MC PAC aIW aW twin;'])

end
%% confusion matrix
%clear
pid='130';
eval(['load D:\Tone\',pid,'\',pid,'_classification_result;']) ;  K=6;
ma=max(MC(:));
for n=1:16
    figure(24); imagesc(MC(:,:,n)); set(gca,'ydir','norm','fontsize',14,'fontweight','bold')
    load cm9; colormap(cm9);caxis([1/K+0.02  ma]);colorbar;
    ylabel(' Presented stim'); xlabel(' Predicted stim');
    set(gcf,'paperunits','inches','papersize',[6 5],...
 'PaperPositionMode','manual','paperunits','normalized','paperposition',[0 0 1 1]);
    eval(['print -f24 -depsc -r150 D:\Tone\',pid,'\',pid,'_cm_tw',num2str(n),';'])
    eval(['print -f24 -dpng -r300 D:\Tone\',pid,'\',pid,'_cm_tw',num2str(n),';'])
end
%%  adjusted weights image
clear
pid='130';
eval(['load D:\Tone\',pid,'\',pid,'_classification_result;']) ;  K=6;
dat=AW.*repmat((Accuracy(1,:)/100-1/K),size(AW,1),1);

    dat1=AW(:,[1:3 15:16]).*repmat((Accuracy(1,[1:3 15:16])/100-1/K),size(AW,1),1);
    cdc=dat1(:);
    ff=find(cdc~=0);
    [a]=prctile(abs(cdc(ff)),[90]);
    dat1=AW(:,[4:14]).*repmat((Accuracy(1,[4:14])/100-1/K),size(AW,1),1);
    cdc=dat1(:);
    ff=find(cdc~=0);
    [b]=prctile(abs(cdc(ff)),[95]);

g=meshgrid(1:12,1:8);g=g(:);  g2=repmat([1:8],1,12) ;
clear cdc dat1 ff 
for n=1:16
    figure(30); imagesc(reshape(dat(:,n),8,12)); set(gca,'ydir','norm','fontsize',8,'fontweight','bold','color',[0.5 0.5 0.5])
    load cm11; colormap(cm11);caxis([a  b]);colorbar;
    title(' Weights mapped on 96-contacts surface grid');
    for ii=1:96
     eval(['h=text(g(ii)-0.15,g2(ii),''',num2str(ii),''');']); set(h,'color',[0.8 0.85 0.85]);
    end
    set(gcf,'paperunits','inches','papersize',[8 4.2],...
 'PaperPositionMode','manual','paperunits','normalized','paperposition',[0 0 1 1]);
    eval(['print -f30 -depsc -r150 D:\Tone\',pid,'\',pid,'_aweights_tw',num2str(n),';'])
    eval(['print -f30 -dpng -r300 D:\Tone\',pid,'\',pid,'_aweights_tw',num2str(n),';'])
end
    twin=[-0.1:0.05:0.65]+0.025;
    
    figure(100);plot(twin,Accuracy(1,:),'.-','linewidth',2);grid on;
    xlabel('sec');ylabel('percent accuracy');
    set(gca,'fontsize',12,'fontweight','bold')
    eval(['print -f100 -dpng -r300 D:\Tone\',pid,'\',pid,'_accuracy;'])
     
eval(['save D:\Tone\',pid,'\',pid,'_classification_result MW AW  IW Accuracy MC PAC aIW aW twin dat;'])