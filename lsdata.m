clear;
clc;
%固定数值
e=1.6*10^-19;
ee=8.85*10^-12;      %电介常数
u=4*pi*10^-7;
M=2*1.67*10^-27;
m=9.1*10^-31;
w=4.6*2*pi*10^9; %低杂波波源频率
c=299792458;


%输入量
load('profiles_080307_004000.mat'); 
rho=[0:0.01:1];
a=0.45;
I=400*10^3;
It=11000;
B0=It*1.7/4085/1.85;
B=B0*1.85./(1.85+a.*rho);
deltan=0.03; 
Npar=2.04;
Z=2;


Te=Te*1000*e;
Ti=Te./TeTi;
ni=ne;




wpe=(ne.*e*e/ee/m).^0.5;   %电子等离子体频率
wce=e.*B/m;  %电子回旋频率
wpi=(ni.*e*e/ee/M).^0.5;
wci=e.*B/M;

kpar=Npar*w/c; %k//值

%kesi0=2*k0*sind(sitas/2 )

aa=wpe./wce;
bb=wpi./w;
npar_acc=wpe./wce+(1+(wpe./wce).^2-(wpi./w).^2).^(0.5);  %可近性条件

%wpi,wpe,wce
w_LH=(wpi.^2./(1+(wpe./wce).^2)).^0.5;  %LH共振频率
 
w0=w./w_LH;
% beta=(ni*Ti+ne*Te)/(B)^2*2*u;
%betan=(beta*100)./(a*B/I/u);    %betan 可直接从实验中得出，可以反推beta大小
%betan=0.458;
%beta=betan*(a*B/I/u)/100;
%krou=Npar*(beta*M/2/m)^0.5*w0/(w0^2-1)^0.5; %berger给出的kperp
roui=(Z/2*M*Ti).^0.5./(e.*B);
%kperp=krou/roui;
%k0=kperp;
k0=sqrt(wpe.^2*kpar^2/w^2);    %基本上无误，可由β反推
ksi0roui=0.09;
ksi01=ksi0roui./roui;
%ksi0=(m/M*m*c^2./Ti*(1+wce^2./wpe.^2)*ksi0roui.^2).^0.5;
ksi0=((1-w0.^-2)./Npar.^2*m/M*m*c^2./Ti.*(1+wce.^2./wpe.^2)*ksi0roui.^2).^0.5.*k0;
ksirat=ksi0./ksi01;
while abs(ksirat)<1.25
    break;
  ksi0roui=ksi0roui+0.01;
end

for rh=[1:1:101]
    err=k0-ksi01(rh);
    if abs(err)<1
        break;
    end
end

%x1=(ksi0(rh)*c./wpe(rh)./Npar).^2*a/2;
x11=B.^2./(2*u*ni(1).*Ti(rh))*m/M.*ksi0roui.^2/Npar.^2*a/2; 
%l=a-x1;
ls=2/pi.*(w0.^2-1).^2./ksi0./w0.^2.*wci.^2./w_LH.^2*deltan.^-2;
%lst=3/16*pi^0.5*deltan.^2*(ksi0./k0).^2*(wpe.^2./w./wce).^2*ksi0;
%lls=l./lst;
%tao=ksi0roui./Npar^2*(wpi(1)./w).^2*(m*c^2/2./Ti).^0.5*c/a./(wce.*wci).^0.5;
tao=1./Npar^2*(wpi(1)./w).^2*(m*c^2/2./Ti).^0.5*c/a./(wce.*wci).^0.5;

