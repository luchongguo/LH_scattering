function [lent,lenl,len,cnzm,cnrm,cmm,beta11,Vgperm,nu1]=Copy_2_of_ks(z,r,phi,cnz,cnr,cm,cnpar,skperp,skperm,ktheta,deltan)

%%
global Bt1 Br1 Bz1 eqdsk_r1 eqdsk_z1 Ni1 Ne1  w Btot1

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

c=299792458;
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

step=10^-5;
hw=step*w;
wp=w+hw;
dwp=w/wp;
xar1p=xar1*dwp*dwp;
xar2p=xar2*dwp*dwp;
yar1p=yar1*dwp;
yar2p=yar2*dwp;
[epsp]=calc_eps1(xar1p,xar2p,yar1p,yar2p);


wm=w-hw;
dwm=w/wm;
xar1m=xar1*dwm*dwm;
xar2m=xar2*dwm*dwm;
yar1m=yar1*dwm;
yar2m=yar2*dwm;
% eps.perpm=1+xar1m/(yar1m)^2-xar2;
% eps.parm=1-xar1m-xar2m;
% eps.xym=xar1m/yar1m;
[epsm]=calc_eps1(xar1m,xar2m,yar1m,yar2m);
% deps.perp=(epsp.perp-epsm.perp)/2/hw;
% deps.par=(epsp.par-epsm.par)/2/hw;
% deps.xy=(epsp.xy-epsm.xy)/2/hw;
% wdepsperp=w*deps.perp;
% wdepspar=w*deps.par;
% wdepsxy=w*deps.xy;
wdepsperp=(wp*epsp.perp-wm*epsm.perp)/2/hw;
wdepspar=(wp*epsp.par-wm*epsm.par)/2/hw;
wdepsxy=(wp*epsp.xy-wm*epsm.xy)/2/hw;
%%      ��H���������������������V
hxx=eps.perp+wdepsperp/2;
hxy=eps.xy+wdepsxy/2;
hzz=eps.par+wdepspar/2;
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

Hh=hxx.*(aa.^2+bb.^2)+2.*aa.*bb.*hxy+hzz;
Hhm=hxx.*(aam.^2+bbm.^2)+2.*aam.*bbm.*hxy+hzz;
%%
e=1.6*10^-19;
ee=8.85*10^-12;      %��鳣��
M=2*1.67*10^-27;
m=9.1*10^-31;
c=299792458;
%w=4.6*2*pi*10^9; %���Ӳ���ԴƵ��

%%
wpe=(ne*e*e/ee/m).^0.5;   %���ӵ�������Ƶ��
wce=e*Btot1(ii,jj)/m;  %���ӻ���Ƶ��
wpi=(ni*e*e/ee/M).^0.5;
wci=e*Btot1(ii,jj)/M;

beta1=-pi:0.01:pi;
%[beta1]=calc_random_beta(uzz,urr,uii,skperp,skperm,ktheta);
for iii=1:length(beta1)


[cnzm,cnrm,cmm]=rotate(cnz,cnr,cm,r,br,bz,bphi,cnpar,skperp,skperm,w,beta1(iii)); % �����µģ�cnz,cnr,cm���µ�Ⱥ�ٶ�
%ddd %�����������ʱ�����
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

Vgperms=(Vgtot2-Vgpar^2);       
Vgperm=sqrt(Vgperms); 
cj11=uzz+cos(beta1(iii)).*urr;
cj12=sin(beta1(iii)).*uii;
cj1=cj11.^2+cj12.^2;

%exp12=exp(-(skperp^2+skperm^2)/ktheta^2);
lenls2=skperp^2+skperm^2-2*skperp*skperm*cos(beta1(iii));
exp3=exp(-lenls2/ktheta^2);
%Vkkm2=0.25*w^2*abs((Hh*Hhm*cj)^-1);     %  |V|����
Vkkm2=0.25*w^2*abs((Hh*Hhm)^-1);     %  |V|����
const=2*pi.*skperm./Vgperm.*Vkkm2.*(pi.*ktheta^2)^-1.*deltan^2;
%pks=const.*cj1.*exp3;  

ktot2=-skperp2*((wpe./wce).^2-(wpi./w).^2)-skpr2.*(-wpe.^2-wpi.^2)./w./w;
% P_ott=(1-eps.perp.*(sin(beta1/2)).^2+(sin(beta1)).^2*wpe.^4/w/w/wce/wce).*...
%         exp(-4*skperp.^2/ktheta.^2.*(sin(beta1/2)).^2);
% P_ott1=(1+2*skperp2./ktot2.*((wpe./wce).^2-(wpi./w).^2).*(sin(beta1/2)).^2+(sin(beta1)).^2*wpe.^4/w/w/wce/wce.*(skperp2./ktot2).^2).*...
%         exp(-4*skperp.^2/ktheta.^2.*(sin(beta1/2)).^2);
P_Bonoli(iii)=const.*((uzz+cos(beta1(iii))*urr).^2+(sin(beta1(iii))*uii).^2).*...
    exp(-(skperp^2+skperm^2-2*skperp*skperm*cos(beta1(iii)))/ktheta^2);
P_Bonoli_1(iii)=P_Bonoli(iii)./const;
nu(iii)=(2*(sin(beta1(iii)/2)).^2*const.*((uzz+cos(beta1(iii))*urr).^2+(sin(beta1(iii))*uii).^2).*...
    exp(-(skperp^2+skperm^2-2*skperp*skperm*cos(beta1(iii)))/ktheta^2));

end
    Ptot_Bonoli=trapz(beta1,P_Bonoli);
    nu1=trapz(beta1,nu);
    
    lent=Ptot_Bonoli;

while 1

    beta2=rand(1)*2*pi-pi;%����[-pi,pi]���ȷֲ������
    [~,ind_beta1]=min(abs(beta2-beta1));
    
    if beta2<=pi
        f=P_Bonoli(ind_beta1)/max(P_Bonoli);

    end         %�����Ӧ�ܶȺ���ֵf(beta2)
    r=rand(1);  %����[0,1]���ȷֲ������
    if r<=f     %��������rС��f(beta2)�����ɸ�t����������a��
        beta11=beta2;
        break;
        
    end
end

%% ����v��
beta11=abs(beta11);
[cnzm,cnrm,cmm]=rotate(cnz,cnr,cm,r,br,bz,bphi,cnpar,skperp,skperm,w,beta11); % �����µģ�cnz,cnr,cm���µ�Ⱥ�ٶ�
[~,ddt]=calc_ddd1_old(z,r,phi,cnzm,cnrm,cmm);   %ddd %�����������ʱ�����
vgz=ddt.z;      %����Ⱥ�ٶȵ���������
vgr=ddt.r;
vgphi=r*ddt.phi;

% Vgdotb=vgz*bz+vgr*br+vgphi*bphi;
btot2=bz^2+br^2+bphi^2;
% Vgperms=(vgz-Vgdotb*bz/btot)^2+(vgr-Vgdotb*br/btot)^2+(vgphi-Vgdotb*bphi/btot)^2
Vgtot2=vgz^2+vgr^2+vgphi^2;
btot=sqrt(btot2);
Vgpar=(vgz*bz+vgr*br+vgphi*bphi)/btot;

Vgperms=(Vgtot2-Vgpar^2);       
Vgperm=sqrt(Vgperms); 
lenl=Vgperm/nu1;
len=0;