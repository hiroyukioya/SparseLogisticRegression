function [cm]=mycolor
%%

cm(:,1) = interp1([1 100 256],[0.1  0.7  1],[1:256],'linear');
cm(:,2) = interp1([1 128 256],[0.1  0.25  1],[1:256],'spline');
cm(:,3) = interp1([1 128 256],[0.1  0.05  0],[1:256],'linear');
plot(cm(:,1),'r');hold on
plot(cm(:,2),'g');hold on
plot(cm(:,3),'b');hold on
save /home/hiroyuki/sfn2010/cm  cm
%%

cm2(:,1) = interp1([1 64 128 192 256],[0   0.3  1  1   1],[1:256],'linear');
cm2(:,2) = interp1([1 64 128 192 256],[0   0.3  1  0.3   0],[1:256],'linear');
cm2(:,3) = interp1([1 64 128 192 256],[1   1  1   0.3  0],[1:256],'linear');
clf;plot(cm2(:,1),'r','linewidth',2);hold on
plot(cm2(:,2),'g','linewidth',2);hold on
plot(cm2(:,3),'b','linewidth',2);hold on
save /home/hiroyuki/sfn2010/cm2  cm2
%%

cm3(:,1) = interp1([1 64 128 192 256],[0   0.25  0.6   0.8     1],[1:256],'linear');
cm3(:,2) = interp1([1 64 128 192 256],[0   0.2    0.6   0.2    0],[1:256],'linear');
cm3(:,3) = interp1([1 64 128 192 256],[1   0.8    0.6   0.25  0],[1:256],'linear');
clf;plot(cm3(:,1),'r','linewidth',2);hold on
plot(cm3(:,2),'g','linewidth',2);hold on
plot(cm3(:,3),'b','linewidth',2);hold on
save /home/hiroyuki/sfn2010/cm3  cm3
%%

cm4(:,1) = interp1([1 40 128 216 256],[0   0.25    1      1     1],[1:256],'linear');
cm4(:,2) = interp1([1 40 128 216 256],[0   0.25    1   0.25    0],[1:256],'linear');
cm4(:,3) = interp1([1 40 128 216 256],[1     1       1   0.25  0],[1:256],'linear');
clf;plot(cm4(:,1),'r','linewidth',2);hold on
plot(cm4(:,2),'g','linewidth',2);hold on
plot(cm4(:,3),'b','linewidth',2);hold on
save /home/hiroyuki/sfn2010/cm4  cm4
%%

cm5(:,1) = interp1([1 40 128 216 245 256],[0   0.1    0.6     0.8    1 1],[1:256],'linear');
cm5(:,2) = interp1([1 40 128 216 245 256],[0    0.1      0.6     0.8   1 1],[1:256],'linear');
cm5(:,3) = interp1([1 40 128 216 245 256],[0    0      0   0     0   0],[1:256],'linear');
save /home/hiroyuki/sfn2010/cm5  cm5;
%%

cm6(:,1) = interp1([1 40 128 216 235 256],[1    1        1            1        1         1],[1:256],'linear');
cm6(:,2) = interp1([1 40 128 216 235 256],[1    0.9      0.5       0.3      0.05      0],[1:256],'linear');
cm6(:,3) = interp1([1 40 128 216 235 256],[1    0.9      0.4       0.1      0.02       0],[1:256],'linear');
save /home/hiroyuki/sfn2010/cm6  cm6;
%% cm7
cm7(:,1) = interp1([1 40 128 216 235 256],[0.05    0.2        0.6           0.7        0.8         1],[1:256],'linear');
cm7(:,2) = interp1([1 40 128 216 235 256],[0.05    0.2          0.6          0.7         0.8       1],[1:256],'linear');
cm7(:,3) = interp1([1 40 128 216 235 256],[0.05    0.05         0.05        0.05       0.05         0],[1:256],'linear');
save /home/hiroyuki/sfn2010/cm7  cm7;
%% cm8
cm8(:,1) = interp1([1 40 128 216 256],[1    0.8        0.8        0.9     1],[1:256],'linear');
cm8(:,2) = interp1([1 40 128 216 256],[1    0.8        0.4       0.25      0.1],[1:256],'linear');
cm8(:,3) = interp1([1 40 128 216 256],[1    0.8        0.3       0.15      0],[1:256],'linear');
save /home/hiroyuki/sfn2010/cm8  cm8
%% cm9
cv=0.95;       val=1;

cm9(:,1) = interp1([1 40 128  200  240 256],[ cv cv cv cv cv cv],[1:256],'linear');
cm9(:,2) = interp1([1 40 128  200  240 256],[0    0.1    0.3    0.8  1   1],[1:256],'linear');
cm9(:,3) = interp1([1 40 128  200  240 256],[val   val    val     val    1 1],[1:256],'linear');

cm9=hsv2rgb(cm9);
save /home/hiroyuki/sfn2010/cm9  cm9

%% cm10
cv=0.155:0.1:0.9;   
cv=linspace(0.05,0.18,6);
val=1;

cm10(:,1) = interp1([1 40 128  200  240 256],[cv(1) cv(2) cv(3) cv(4) cv(5) cv(6)],[1:256],'linear');
cm10(:,2) = interp1([1 40 128  200  240 256],[0.5  0.7   0.8   0.9    1     1],[1:256],'linear');
cm10(:,3) = interp1([1 40 128  200  240 256],[0     0.3   0.65   0.8  0.9  1],[1:256],'linear');

cm10=hsv2rgb(cm10);
%save /home/hiroyuki/sfn2010/cm10  cm10
save D:\SLCtool\cm10  cm10
%%
%% cm11
cv=0.02;
%cv=linspace(0.05,0.18,6);
val=1;

cm11(:,1) = interp1([1 40 128  200  240 256],[cv(1) cv(1) cv(1) cv(1) cv(1) cv(1)],[1:256],'linear');
cm11(:,2) = interp1([1 40 128  200  240 256],[0  0.20   0.5   0.75    0.9     1],[1:256],'cubic');
cm11(:,3) = interp1([1 40 128  200  240 256],[1    1     1     1     0.95    0.9],[1:256],'cubic');

cm11=hsv2rgb(cm11);
%save /home/hiroyuki/sfn2010/cm10  cm10
save D:\SLCtool\cm11  cm11