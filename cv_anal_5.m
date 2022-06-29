%% Load Data
clear;
for patient=1
    if patient==1
        % PT180
        clear sid stime data cutfp
        load ~/matdata/180-066cvL;
        temp=cutfp; f=find(sid~=255); temp=temp(:,f,:); clear cutfp; 
        load ~/matdata/180-066cvH;
        data=cat(3,temp,cutfp(:,f,:));
        sid1=sid(f); stime1=stime(f);
        
        clear temp f cutfp sid ;
        load ~/matdata/180-067cvL;
        temp=cutfp; f=find(sid~=255); temp=temp(:,f,:); clear cutfp; 
        load ~/matdata/180-067cvH;
        data2=cat(3,temp,cutfp(:,f,:));
        sid2=sid(f); stime2=stime(f);
        
        data=cat(2,data,data2);
        sid=[sid1 sid2];   stime=[stime1 stime2]; 
        clear  sid1 sid2 temp data2 f stime1 stime2 cutfp Ntrial T J 
        
        eval(['save ~/cv/',num2str(patient),'/Adata  data LCH HCH  sid stime tstamp ;'])
        clear cutfp Ntrial fs temp f        
    elseif patient==2
        % PT186-1
        clear sid stime data cutfp
        load ~/matdata/186-062cvL;
        temp=cutfp; f=find(sid~=255); temp=temp(:,f,:); clear cutfp; 
        load ~/matdata/186-062cvH;
        data=cat(3,temp,cutfp(:,f,:));
        sid1=sid(f); stime1=stime(f);  
        
        clear temp f cutfp sid ;
        load ~/matdata/186-064cvL;
        temp=cutfp; f=find(sid~=255); temp=temp(:,f,:); clear cutfp; 
        load ~/matdata/186-064cvH;
        data2=cat(3,temp,cutfp(:,f,:));
        sid2=sid(f); stime2=stime(f);
        clear sid stime cutfp temp
        data=cat(2,data,data2);
        sid=[sid1 sid2];   stime=[stime1 stime2]; 
        clear  sid1 sid2 temp data2 f stime1 stime2 cutfp Ntrial T J 
        
        eval(['save ~/cv/',num2str(patient),'/Adata  data LCH HCH  sid stime tstamp;'])
        clear cutfp Ntrial fs temp f              
        
    elseif patient==3
        % PT186-2
        clear sid stime data cutfp
        load ~/matdata/186-067cvL;
        temp=cutfp; f=find(sid<254); temp=temp(:,f,:); clear cutfp; 
        load ~/matdata/186-067cvH;
        data=cat(3,temp,cutfp(:,f,:));
        sid1=sid(f); stime1=stime(f);
        
        clear temp f cutfp sid ;
        load ~/matdata/186-069cvL;
        temp=cutfp; f=find(sid<254); temp=temp(:,f,:); clear cutfp; 
        load ~/matdata/186-069cvH;
        data2=cat(3,temp,cutfp(:,f,:));
        sid2=sid(f); stime2=stime(f);
        
        data=cat(2,data,data2);
        sid=[sid1 sid2];   stime=[stime1 stime2]; 
        clear  sid1 sid2 temp data2 f stime1 stime2 cutfp Ntrial T J 
        
        eval(['save ~/cv/',num2str(patient),'/Adata  data LCH HCH  sid stime tstamp;'])
        clear cutfp Ntrial fs temp f              
        
    elseif patient==4
        % PT198
        clear sid stime data cutfp
        load ~/matdata/198-058cvL;
        temp=cutfp; f=find(sid<254); temp=temp(:,f,:); clear cutfp; 
        load ~/matdata/198-058cvH;
        data=cat(3,temp,cutfp(:,f,:));
        sid1=sid(f); stime1=stime(f);

        clear temp f cutfp sid ;
        load ~/matdata/198-059cvL;
        temp=cutfp; f=find(sid<254); temp=temp(:,f,:); clear cutfp; 
        load ~/matdata/198-059cvH;
        data2=cat(3,temp,cutfp(:,f,:));
        sid2=sid(f); stime2=stime(f);
   
        data=cat(2,data,data2);
        sid=[sid1 sid2];   stime=[stime1 stime2]; 
        clear  sid1 sid2 temp data2 f stime1 stime2 cutfp Ntrial T J 
        
        eval(['save ~/cv/',num2str(patient),'/Adata  data LCH HCH  sid stime tstamp;'])
        clear cutfp Ntrial fs temp f              

    end
%%   Do TMT.........
    st=find(tstamp<=2);
    data=data(st,:,:);
    tstamp=tstamp(st);
    clear st
    cd ~/cv
    %   set-up parameters .................
    inparg.FS=FS;   inparg.taper=[2 3];
    inparg.wlength=128;  inparg.stepsize=25;
    inparg.preT=1;   inparg.reftime=[-0.85  -0.25];  inparg.sndtrial=[];
    inparg.reftrial=[]; inparg.mode=1; inparg.figrange=[-0.8 1.8];
    % Do MTM spectral analysis .............
    for ch=1:size(data,3)
        close all;
        [Z, S , TindH,  fax,  FX,T,J]=TMTspectrogram(data(:,:,ch), inparg);   
        set(gcf,'paperunits','inches','papersize',[7 4],...
     'PaperPositionMode','manual','paperunits','normalized','paperposition',[0 0 1 1]);
        eval(['print -f1 -dpng -r300 ~/cv/',num2str(patient),'/tmt-',num2str(ch),'.png;']); 
        eval(['save ~/cv/',num2str(patient),'/fdata-',num2str(ch),'  sid stime Z S TindH T J fax ;'])
    end
    clear T J gdata FX;
end
%%
%  /////////////**********CLASSIFICATION***********///////////////
%  //////////////   Make input feature matrix     ////////////////////
%  ///////////////   Z = Time x trial x ch x band    /////////////////

for patient=1
    cd /home/hiroyuki/cv
    clear sid x Z V
    eval(['load ~/cv/',num2str(patient),'/Adata;'])
    eval(['load ~/cv/',num2str(patient),'/fdata-',num2str(ch),';']); 
%     passband=[2 12;  15 25 ; 70 140];
    passband=[70 150];
    st=find(tstamp<=2);
    data=data(st,:,:);
    tstamp=tstamp(st);
    totalch=size(data,3);
    Tind=[-1:0.05:2];
    Z=zeros(length(Tind)-1, size(data,2), totalch,1);
    Rejm=zeros( size(data,2),totalch,3);
    clear st stime exp 
    
    %---------------------------------------
    for t=1:length(Tind)-1
         tt=find(tstamp>Tind(t)-0.015 & tstamp<=Tind(t)+0.015);
         V(t,:,:)=squeeze(mean(data(tt,:,:),1));
    end
    
    % Filtering  --------------------------------
    Zdat=[]; 
    for band=1:size(passband,1)
        for ch=1:size(data,3);
             eval(['load ~/cv/',num2str(patient),'/fdata-',num2str(ch),';']);clear FX S gdata
             f2=find(fax>=70 & fax<=140); 
             p2=squeeze(mean(Z(f2,:,:),1)); 
             Rejm(J,ch)=1;
             Zdat(:,:,ch)=p2';
        end
    end
    clear band ch t tt Zt 
    
    SEQ=sid;
    rejtrial=find(sum(Rejm(:,:,1),2)>size(Z,1)*0.2);   % Reject bad trials..........
    SEQ(rejtrial)=[];    Z(:,rejtrial,:,:)=[]; V(:,rejtrial,:)=[];
    fprintf(' Number of discarded trials = %2.0f\n',length(rejtrial));
    
    % Prefiltering ...........................................................
    mv=squeeze(mean(V,2)); mmv=max(abs(mv),[],1);
    u=find(Tind>0 & Tind<2);
    mx=squeeze(mean(Z,2)); mxmi=squeeze(min(mx,[],1));  mxma=squeeze(max(mx,[],1));
    if size(passband,1)==1
        f=find(mxma>4);
        f2=find(mxmi<-4);
        f3=find(mmv>40);
        selchhg=union(f,f2);
        selchlw=f3+totalch;
        selc=[selchhg selchlw];
    else
        for n=1:size(passband,1)
            f=find(mxma(:,n)>4);
            f2=find(mxmi(:,n)<-4);
            eval(['ff',num2str(n),'=union(f,f2);'])
        end 
        selc=[ff1 ;ff1+totalch ;ff2+totalch*2];
    end
    fprintf(' Number of selected channels = %2.0f\n',length(selc));
    clear ff* f f2 mx* u prestim rejp* a1 a2 ans N J T Fc* 

    %%  /**********   Classification    ***********/
    u=unique(SEQ);
    sortseq=[];
     for n=1:length(u)
         o=find(SEQ==u(n));
         sortseq(o)=n;
     end
    [JT,ii]=sort(sortseq);
    clear  n o sortseq  o n sortseq lch hch tch J 
        
    %  Do classification ================
    cd ~/SLCtool; 
    
    matlabpool open 2
    adW=[]; MW=[]; MC=[] ; PAC=[]; accuracy=[];AC=[];
    for n=20:55
%          z=[squeeze(Z(n,:,:,1)) squeeze(Z(n,:,:,2)) squeeze(Z(n,:,:,3))];z=z(:,selc);
         z=[squeeze(Z(n,:,:)) squeeze(V(n,:,:))];z=z(:,selc);
         datt=centernormalize(z,1);
         dat=datt(ii,:);
         [adW(:,n), MW(:,n), W, AC, MC(:,:,n), PAC(:,:,n)] = SLCtool(dat, (JT)');
         accuracy(n)=AC.ec; accuracy2(n)=AC.pc;
         tw=Tind(n),
         clear dat
         eval(['print -f24 -dpng -r300 ~/cv/',num2str(patient),'/CM',num2str(n),';'])
         eval(['print -f26 -dpng -r300 ~/cv/',num2str(patient),'/PCM',num2str(n),';'])
         eval(['saveas(figure(24),''~/cv/',num2str(patient),'/CM',num2str(n),''');'])
         eval(['saveas(figure(26),''~/cv/',num2str(patient),'/PCM',num2str(n),''');'])
    end
  matlabpool close;
  
  adW(1,:)=[]; MW(1,:)=[];
  Lweights=zeros(totalch,size(adW,2));
  Hweights=zeros(totalch,size(adW,2));
  
  gg=find(selc<=totalch);
  Zweights(selc(gg),:)=adW(gg,:);
  gg2=find(selc>totalch);
  Vweights(selc(gg2)-totalch,:)=adW(gg2,:);
  
  eval(['save ~/cv/',num2str(patient),'/slcdata  Tind selc adW MW W AC MC PAC accu* z Hweights Lweights  ;'])
end
%%

% Combine images into one big matrix
close all;
C=[8:8:96 7:8:96 6:8:96 5:8:96 4:8:96 3:8:96 2:8:96 1:8:96];
X=[];X1=[];X2=[];X3=[];X4=[];X5=[];X6=[];X7=[]; X8=[];
for n=1:96
A=[];
eval(['A=imread(''~/cv/1/tmt-',num2str(C(n)),'.png'');'])
A(: ,[1:309 1756:end] ,:)=[]; A([1:123 624:end],:,:)=[];
A=imresize(A,0.5);
[a,b,c]=size(A);
A=cat(2,A,uint8(zeros(a,10,c)));
    if n<=12
    X1=[X1 A];
    elseif n>=13 & n<=24
    X2=[X2 A];
    elseif n>=25 & n<=36
    X3=[X3 A];
    elseif n>=37 & n<=48
    X4=[X4 A];
    elseif n>=49 & n<=60
    X5=[X5 A];
    elseif n>=61 & n<=72
    X6=[X6 A];
    elseif n>=73 & n<=84
    X7=[X7 A];
    elseif n>=85 & n<=96
    X8=[X8 A];
    end
end
X=[X1;X2;X3;X4;X5;X6;X7;X8]; clear X1 X2 X3 X4 X5 X6
imshow(X);
set(gcf,'paperunits','inches','papersize',[12 6],'PaperOrientation', 'portrait',...
'PaperPositionMode','manual','paperunits','normalized','paperposition',[0 0 1 1]);
eval(['print -f1 -dpng -r600 ~/cv/1/TMT-combined_LTW;'])

temp=Hweights(1:96,24)';
imagesc(flipud(reshape(temp,8,12)));

temp=Lweights(1:96,24)';
imagesc(flipud(reshape(temp,8,12)));
