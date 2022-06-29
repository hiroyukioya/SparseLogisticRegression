clear;

pid=178
if pid==166
    block='017'; blk=1; reso=0;
    aa=[0 4 16 28 40 44];Nn=3;
    NF=[aa aa+1 aa+2; 1:6 1:6 1:6];
elseif  pid==175
    block='023'; blk=1; reso=0;
    aa=[0 6 24 42 60 66];Nn=4;
    NF=[aa aa+1 aa+2  aa+3; 1:6 1:6 1:6 1:6];
elseif  pid==178
    block='018'; blk=1; reso=0;
    aa=[0 6 24 42 60 66];Nn=4;
    NF=[aa aa+1 aa+2  aa+3; 1:6 1:6 1:6 1:6];
end
pid=num2str(pid);
    
for w=1:16
      ww=w;
      eval(['nx',num2str(ww),'=[];']);eval(['nnx',num2str(ww),'=[];']);eval(['Rnx',num2str(ww),'=[];']);
      T=[];ntr=0;
        for n=1:Nn
           eval(['nx',num2str(ww),'=[];']);
           for fr=(n-1)*6+1:n*6
                stimfreq=NF(2,fr);
                stimid=NF(1,fr);
                eval(['a1=importdata(''D:\Tone\GT8s\All_trials_ERBP_',pid,'-',block,'_blk',num2str(blk),'_',num2str(stimid),'_trials_win',num2str(w),'.txt'');'])
                sx1x=a1(:,[2:end])';    ch=a1(:,1); 
                clear a1; [a,b]=size(sx1x);
                eval(['nx',num2str(ww),'=[nx',num2str(ww),';sx1x];']);   clear sx1x;
                tempT=ones(a,1)*stimfreq;  T=[T; tempT];
           end     
            eval(['nnx',num2str(ww),'=centernormalize(nx',num2str(ww),',1);']);
            eval(['Rnx',num2str(ww),'=[Rnx',num2str(ww),'; nnx',num2str(ww),'];']);
            clear nx*
        end
       eval(['x',num2str(ww),'=centernormalize(Rnx',num2str(ww),',1);']);
end
 
    twin=[-0.1:0.05:0.65]+0.025;
    
    
 %.......................... for saving nx* (non-normalized data) .......................
 for w=1:16
      ww=w;
      eval(['nx',num2str(ww),'=[];']);
      eval(['nx',num2str(ww),'=[];']);
      for n=1:Nn
           for fr=(n-1)*6+1:n*6
                stimfreq=NF(2,fr);
                stimid=NF(1,fr);
                eval(['a1=importdata(''D:\Tone\GT8s\All_trials_ERBP_',pid,'-',block,'_blk',num2str(blk),'_',num2str(stimid),'_trials_win',num2str(w),'.txt'');'])
                sx1x=a1(:,[2:end])';    ch=a1(:,1); 
                clear a1; [a,b]=size(sx1x);
                eval(['nx',num2str(ww),'=[nx',num2str(ww),';sx1x];']);   clear sx1x;
           end    
      end
 end   
 
 % Sort data -------
 [T,I]=sort(T);
 for w=1:16
     eval(['x',num2str(w),'=x',num2str(w),'(I,:);'])
 end
     
clear  a temp i n s s2 tempT b n w ww i temp* x nxn*  nnx* a1 n fr Rnx*
eval(['save D:\Tone\GridToneAnal_pt',pid,'  T x* ch nx* twin;']);
