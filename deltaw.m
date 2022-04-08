%固定数值
e=1.6*10^-19;
ee=8.85*10^-12;      %电介常数
M=2*1.67*10^-27;
m=9.1*10^-31;
w=4.6*2*pi*10^9; %低杂波波源频率
c=299792458;
a=0.45;

Npar=4;
B=6;
ne=2e20;
ni=2e18;
Ti=50*1.6e-19;
wpe=(ne.*e*e/ee/m).^0.5;   %电子等离子体频率
wce=e*B/m;  %电子回旋频率
kpar=Npar*w/c;
kperp=wpe/w*kpar;
wpi=(ni.*e*e/ee/M).^0.5;
wci=e*B/M;
roui=(Ti*2/M)^0.5/wci;
deltaw1=Npar/2*(c/a)*(wpe/wce)*(2*Ti/m/c^2);
deltaw2=(pi)^0.5*(c/a)^2/6*(wpi/w)^2*(wpe/wce)^2*(2*Ti/m/c^2)^1.5*(m/M)^0.5*(c/a/wci);
detlanpar=0.5*(1+wpe^2/wce^2)^0.5*(M*c*c/Ti)*(wci/wpe)^2*0.2*0.2;
detlanpar1=0.5*(1+wpe^2/wce^2)^0.5*(c*0.15/wpe)^2;