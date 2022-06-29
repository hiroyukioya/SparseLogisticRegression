
pid='180';
basedir=strcat('~/Tone/HGtone/', pid)
intens=1;
for ch=1:14
    eval(['load ',basedir,'/SSI_',pid,'_fb1_ch',num2str(ch),'_int',num2str(intens),';'])