function [lent,lenl,len,cnzm,cnrm,cmm,beta1,Vgperm,nu1]=Copy_of_ks(z,r,phi,cnz,cnr,cm,cnpar,skperp,skperm,ktheta,deltan,kmean)

%%
global Bt1 Br1 Bz1 eqdsk_r1 eqdsk_z1 Ni1 Ne1 dNe1 dBrr1 dBzz1 w

% z=ray.z;
% r=ray.r;
% phi=ray.phi;
% cnz=ray.cnz;
% cnr=ray.cnr;
% cm=ray.cm;
[~,ii]=min(abs(z-eqdsk_z1));
[~,jj]=min(abs(r-eqdsk_r1));

br=Br1(ii,jj);
bz=Bz1(ii,jj);
bphi=Bt1(ii,jj);
ne=Ne1(ii,jj);
ni=Ne1(ii,jj);
% dne=dNe1(ii,jj);
% dbrr=dBrr1(ii,jj);
% dbzz=dBzz1(ii,jj);
%%
btot=sqrt(br.^2+bz.^2+bphi.^2);

e=1.6*10^-19;
ee=8.85*10^-12;      %电介常数
M=2*1.67*10^-27;
m=9.1*10^-31;
c=299792458;
wpe=(ne*e*e/ee/m).^0.5;   %电子等离子体频率
wce=e*btot/m;  %电子回旋频率
wpi=(ni*e*e/ee/M).^0.5;

% skperp=snperp*w/c;
% skperm=snperm*w/c;
skperp2=skperp^2;
skperm2=skperm^2;


btot=sqrt(br^2+bz^2+bphi^2);
%skpr=(cnz*bz+cnr*br+cm/r*bphi)/btot;
skpr=cnpar*w/c;

%snpr=(cnz*bz+cnr*br+cm/r*bphi)/btot;
skpr2=skpr^2;

% exx=eps.perp;
% exy=eps,xy;
% ezz=eps.par;
%cosbate=skperp2+skperm2-skperp*skperm


%%      ��wdepsdw����wdD/dw;
[xar1,xar2,yar1,yar2]=calc_xyar(ne,ni,btot);
% eps.perp=1+xar1/(yar1)^2-xar2;
% eps.par=1-xar1-xar2;
% eps.xy=xar1/yar1;
[eps]=calc_eps1(xar1,xar2,yar1,yar2);

uxx=eps.perp-1;
uxy=eps.xy;
uzz=eps.par-1;

hxx=1+wpe.^2./wce.^2;
hxy=0.5*wpe.^2./(w.*wce);
hzz=1;
% step=10^-5;
% hw=step*w;
% wp=w+hw;
% dwp=w/wp;
% xar1p=xar1*dwp*dwp;
% xar2p=xar2*dwp*dwp;
% yar1p=yar1*dwp;
% yar2p=yar2*dwp;
% [epsp]=calc_eps1(xar1p,xar2p,yar1p,yar2p);
% 
% 
% wm=w-hw;
% dwm=w/wm;
% xar1m=xar1*dwm*dwm;
% xar2m=xar2*dwm*dwm;
% yar1m=yar1*dwm;
% yar2m=yar2*dwm;
% % eps.perpm=1+xar1m/(yar1m)^2-xar2;
% % eps.parm=1-xar1m-xar2m;
% % eps.xym=xar1m/yar1m;
% [epsm]=calc_eps1(xar1m,xar2m,yar1m,yar2m);
% % deps.perp=(epsp.perp-epsm.perp)/2/hw;
% % deps.par=(epsp.par-epsm.par)/2/hw;
% % deps.xy=(epsp.xy-epsm.xy)/2/hw;
% % wdepsperp=w*deps.perp;
% % wdepspar=w*deps.par;
% % wdepsxy=w*deps.xy;
% wdepsperp=(wp*epsp.perp-wm*epsm.perp)/2/hw;
% wdepspar=(wp*epsp.par-wm*epsm.par)/2/hw;
% wdepsxy=(wp*epsp.xy-wm*epsm.xy)/2/hw;
% %%      ��H���������������������V
% hxx=eps.perp+wdepsperp/2;
% hxy=eps.xy+wdepsxy/2;
% hzz=eps.par+wdepspar/2;
%%
dtor=(skpr2-eps.perp)*(skpr2+skperp2-eps.perp)-eps.xy^2;
dtorm=(skpr2-eps.perp)*(skpr2+skperm2-eps.perp)-eps.xy^2;

aa=skpr*skperp*(skpr2+skperp2-eps.perp)/dtor;
aam=skpr*skperm*(skpr2+skperm2-eps.perp)/dtorm;
bb=skpr*skperp*eps.xy/dtor;
bbm=skpr*skperm*eps.xy/dtorm;

u1=aa*aam+bb*bbm;
u2=aa*bbm+aam*bb;
urr=uxx*u1+uxy*u2;
uii=uxy*u1+uxx*u2;
%%

%[beta1]=calc_random_beta(uzz,urr,uii,skperp,skperm,ktheta);
[beta1]=Copy_of_calc_random_beta(uzz,urr,uii,skperp,skperm,ktheta,kmean);
%  cosbeta=(skperp^2+skperm^2-ktheta^2)./(2*skperp*skperm);
%  beta222=acos(cosbeta);
%  beta1=beta1+beta222;
%beta1=abs(beta1);
cj11=uzz+cos(beta1)*urr;
cj12=sin(beta1)*uii;
cj1=cj11^2+cj12^2;


Hh=hxx*(aa^2+bb^2)+2*aa*bb*hxy+hzz;
Hhm=hxx*(aam^2+bbm^2)+2*aam*bbm*hxy+hzz;



[cnzm,cnrm,cmm]=rotate(cnz,cnr,cm,r,br,bz,bphi,cnpar,skperp,skperm,w,beta1); % �����µģ�cnz,cnr,cm���µ�Ⱥ�ٶ�
%ddd %�����������ʱ�����
cnzm=real(cnzm);
cnrm=real(cnrm);
cmm=real(cmm);
[~,ddt]=calc_ddd1_old(z,r,phi,cnzm,cnrm,cmm);

vgz=ddt.z;      %����Ⱥ�ٶȵ���������
vgr=ddt.r;
vgphi=r*ddt.phi;

% Vgdotb=vgz*bz+vgr*br+vgphi*bphi;
 btot2=bz^2+br^2+bphi^2;
% Vgperms=(vgz-Vgdotb*bz/btot)^2+(vgr-Vgdotb*br/btot)^2+(vgphi-Vgdotb*bphi/btot)^2
Vgtot2=vgz^2+vgr^2+vgphi^2;
btot=sqrt(btot2);
Vgpar=(vgz*bz+vgr*br+vgphi*bphi)/btot;

Vgperms=abs((Vgtot2-Vgpar^2));       
Vgperm=sqrt(Vgperms);           


%exp12=exp(-(skperp^2+skperm^2)/ktheta^2);
lenls2=skperp^2+skperm^2-2*skperp*skperm*cos(beta1);
exp3=exp(-lenls2/ktheta^2);
%Vkkm2=0.25*w^2*abs((Hh*Hhm*cj)^-1);     %  |V|����
Vkkm2=0.25*w^2*abs((Hh*Hhm)^-1);     %  |V|����
const=abs(2*pi*skperm/Vgperm*Vkkm2*(pi*ktheta^2)^-1*deltan^2);
pks=const.*cj1*exp3;  


beta22=-pi:0.01:pi;
kmean2=kmean*200*pi;


lenls=sqrt(skperp^2+skperm^2-2*skperp*skperm*cos(beta22));
for iii=1:length(beta22)
    if beta22(iii)<0
        lenls(iii)=-abs(lenls(iii));
    else
        lenls(iii)=abs(lenls(iii));
    end
end
pkst1=const.*((uzz+cos(beta22)*urr).^2+(sin(beta22)*uii).^2).*exp(-(lenls-kmean).^2/ktheta^2);
pkstot=trapz(beta22,pkst1);

%pkst1=@(beta2)(const./cj2*exp(-(skperp^2+skperm^2-2*skperp*skperm*cos(beta2))/ktheta^2)); %P��beta��
%pkst1=@(beta2)(const.*((uzz+cos(beta2)*urr).^2+(sin(beta2)*uii).^2).*exp(-(skperp^2+skperm^2-2*skperp*skperm*cos(beta2))/ktheta^2)); %P��beta��
%pkst=const.*cj1.*exp(-(skperp^2+skperm^2-2*skperp*skperm*cos(beta1))/ktheta^2);
%pkstot=integral(pkst1,-pi,pi);        %P��beta����beta�ڣ�-pi��pi���ϵĻ���

if isnan(pkstot) || imag(pkstot)~=0
    pkstot=0;
end

nu=@(beta3)(2*(sin(beta3/2)).^2*const.*((uzz+cos(beta3)*urr).^2+(sin(beta3)*uii).^2).*exp(-(skperp^2+skperm^2-2*skperp*skperm*cos(beta3))/ktheta^2));
nu1=integral(nu,-pi,pi);
if isnan(nu1) || imag(nu1)~=0
    nu1=0;
    lenl=0;    % ls��ֵ
else
    lenl=Vgperm/nu1;    % ls��ֵ
end

lent=pkstot;        % P��beta����beta�ڣ�-pi��pi���ϵĻ���
len=pks;           % K*S