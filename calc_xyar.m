function [xar1,xar2,yar1,yar2,npar_acc]=calc_xyar(ne,ni,B)
%%
global w
%固定数值
e=1.6*10^-19;
ee=8.85*10^-12;      %电介常数
M=2*1.67*10^-27;
m=9.1*10^-31;
c=299792458;
%w=4.6*2*pi*10^9; %低杂波波源频率

%%
wpe=(ne*e*e/ee/m).^0.5;   %电子等离子体频率
wce=e*B/m;  %电子回旋频率
wpi=(ni*e*e/ee/M).^0.5;
wci=e*B/M;
%%
npar_acc=wpe./wce+(1+(wpe./wce).^2-(wpi./w).^2).^(0.5); %可近性条件
xar1=(wpe./w).^2;
yar1=wce./w;
xar2=(wpi./w).^2;
yar2=wci./w;


