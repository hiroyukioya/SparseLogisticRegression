%%
clear;
for patient=4
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
        
        eval(['save ~/cv/',num2str(patient),'/seq sid stime ;'])
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
        
        eval(['save ~/cv/',num2str(patient),'/seq sid stime ;'])
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
        
        eval(['save ~/cv/',num2str(patient),'/seq sid stime ;'])
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
        
        eval(['save ~/cv/',num2str(patient),'/seq sid stime ;'])
        clear cutfp Ntrial fs temp f              

    end

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

%%  /////////////**********CLASSIFICATION***********///////////////
%  //////////////   Make input feature matrix     ////////////////////
%  ///////////////   Z = freqax x trial x TimeW    ///////////////////
clear ;
for patient=3
    clear sid x SEQ 
    eval(['mkdir ~/cv/',num2str(patient),'/RKLC;'])
    totalch=114;
    ch=9;
    eval(['load ~/cv/',num2str(patient),'/fdata-',num2str(ch),';']); 
%     eval(['load ~/cv/',num2str(patient),'/button;']); sid=button; % Subject's decision 
        eval(['load ~/cv/',num2str(patient),'/seq;'])
    SEQ=sid;  
    Rejm=zeros(totalch,size(Z,2));
    totaltrial=size(Z,2);
    clear FX S Z gdata x sid
   
    x=zeros(totaltrial,length(TindH),totalch*2);
    for ch=1:totalch
        pp=[];ch;
        eval(['load ~/cv/',num2str(patient),'/fdata-',num2str(ch),';']);clear FX S gdata
        f1=find(fax>=8 & fax<=20); f2=find(fax>=70 & fax<=140); 
        p1=squeeze(mean(Z(f1,:,:),1)); p2=squeeze(mean(Z(f2,:,:),1));
        pp=cat(3,p1,p2);
        x(:,:,(ch-1)*2+1:ch*2)=pp;
        Rejm(ch,J)=1;  % Rejection matrix
    end
    clear FX S f1 f2 p1 p2 pp ch ans data gd* Z f atem
    x=permute(x,[1 3 2]);
    x=x(:,[1:2:size(x,2) 2:2:size(x,2)],:);  % INPUT MATRIX : x = [ Trial   x   Features  x  TimeWin] 
    if patient==4
%     x([173 427],:,:)=[];       %<---- Rejection----->
    clear sid
    end
    % Prefiltering ...........................................
    u=find(TindH>0 & TindH<2);
    temp1=squeeze(mean(x(:,:,u),1));   
    a=max(temp1(1:totalch,:),[],2);           thre1=max(max(a)/14,2);
    b=min(temp1(1:totalch,:),[],2);            thre2=min(min(b)/14,-2);
    a=max(temp1(totalch+1:end,:),[],2);    thre3=max(max(a)/14,2);
    b=min(temp1(totalch+1:end,:),[],2);     thre4=min(b/4,1);
    clear a*   
    % Determine thresholds .................................
    a1=find(max(temp1(1:totalch,:),[],2)>=thre1 );
    a2=find(min(temp1(1:totalch,:),[],2)<=thre2);
    a3=find(max(temp1(totalch+1:end,:),[],2)>=thre3 ); a3=a3+totalch;
    a4=find(min(temp1(totalch+1:end,:),[],2)<=thre4);   a4=a4+totalch;
    
    selc=union(union(a1,a2),a3);  % No negative gamma band .......

    if patient>=2
        r=[101:107 215:221];    % rejects HG microcontacts outside of HG
        selc=setdiff(selc,r);
    end
    jem=sum(Rejm,2);
    rr=find(jem>totaltrial*0.1);  % Reject bad channels ..............
    selc=setdiff(selc,rr);
    tem=find(selc>totalch); hch=selc(tem)-totalch;
    tem=find(selc<=totalch); lch=selc(tem); tch=union(lch,hch);
    ss=sum(Rejm(tch,:),1);
    sss=find(ss>=5);
    % Reject bad trials ...........................................
    SEQ(sss)=[];
    x(sss,:,:)=[];
    fprintf(' Number of selected channels = %2.0f\n',length(selc));
    fprintf(' Number of discarded trials = %2.0f\n',length(sss));
    clear r a b a1 a2 thre* a*  jem rr  a1 a2 u tempV tch ss sss tem rr jem ii i temp1...
        sid stime u

    
    %  Sort sequence ...........................................
    u=unique(SEQ);
     sortseq=[];
     for n=1:length(u)
         o=find(SEQ==u(n));
         sortseq(o)=n;
     end
    [JT,ii]=sort(sortseq);
    clear  n o sortseq  o n sortseq lch hch tch J 
        
    %  Do classification ================
    cd ~/SLCtool; %matlabpool open 2
    LL=length(u);
    adW=[]; MW=[]; MC=[] ; PAC=[]; accuracy=[];AC=[];
    for n=10:45
         z=x(:,selc,n);
         datt=centernormalize(z,1);
         dat=datt(ii,:);
%          [adW(:,n), MW(:,n), W, AC, MC(:,:,n), PAC(:,:,n)] = SLCtool(dat, (JT)');
         [adW, MW, W, AC, MC, PAC] = RKLCtool(dat, JT', 'linear', 0);

         accuracy(n)=AC.ec; accuracy2(n)=AC.pc;
         tw=TindH(n);
         clear dat
         if LL~=2
             eval(['print -f24 -dpng -r300 ~/cv/',num2str(patient),'/RKLC/CM',num2str(n),';'])
             eval(['print -f26 -dpng -r300 ~/cv/',num2str(patient),'/RKLC/PCM',num2str(n),';'])
             eval(['saveas(figure(24),''~/cv/',num2str(patient),'/RKLC/CM',num2str(n),''');'])
             eval(['saveas(figure(26),''~/cv/',num2str(patient),'/RKLC/PCM',num2str(n),''');'])
         end
    end
  % matlabpool close;
  adW(1,:)=[]; MW(1,:)=[];
  Lweights=zeros(totalch,size(adW,2));
  Hweights=zeros(totalch,size(adW,2));
  gg=find(selc<=totalch);
  Lweights(selc(gg),:)=adW(gg,:);
  gg2=find(selc>totalch);
  Hweights(selc(gg2)-totalch,:)=adW(gg2,:);
  
  eval(['save ~/cv/',num2str(patient),'/RKLC/slcdata  TindH selc adW MW W AC MC PAC accu* x Hweights Lweights  ;'])
end