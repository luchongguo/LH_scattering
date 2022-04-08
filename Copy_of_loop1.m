clear;
clc;
%%
%固定数值

e=1.6*10^-19;
ee=8.85*10^-12;      %电介常数
M=2*1.67*10^-27;
m=9.1*10^-31;
w1=4.6*2*pi*10^9; %低杂波波源频率
c=299792458;

%%
%输入量 b
% 
tstep=1e-11;
nray=4;
lengthray=100;
%%
testDir='F:\cbwu\办公\density_fluction\Gfile';
%gfile = fullfile(testDir,'g080307.040000_old');
gfile = fullfile(testDir,'g080307.040000_old');
gfileVar = readGfile(gfile); 
gvar = gfileVar;
X = linspace(0, gvar.rdim, gvar.nw);
Z = linspace(0, gvar.zdim, gvar.nw);
rgrid = X+gvar.rleft;
zgrid = Z-(gvar.zmid+gvar.zdim/2);
figure(1);
contourf(rgrid,zgrid,-gvar.psirz',50);
hold on;
plot(gvar.rbbbs,gvar.zbbbs,'r','LineWidth',2);   %最外闭合磁面
plot(gvar.rlim,gvar.zlim,'b','LineWidth',2);    %限制器
colorbar();
axis equal;
xlabel('R(m)');
ylabel('Z(m)');
%%
load('nedr.mat');
%load('zrphi.mat');
%load('S48888B.mat');
load('F:\cbwu\办公\density_fluction\B80307\B080307.mat');
%load('ne3.mat');
load('profiles_080307_004000.mat');
[dpsidr,dpsidz]=gradient(gvar.psirz');
X = linspace(0, gvar.rdim, gvar.nw);
Z = linspace(0, gvar.zdim, gvar.nw);
eqdsk_r = X+gvar.rleft;
eqdsk_z = Z-(gvar.zmid+gvar.zdim/2);
[eqdsk_r1]=calc_interp(eqdsk_r);
[eqdsk_z1]=calc_interp(eqdsk_z);
rhopsi=(gvar.psirz'-gvar.simag)/(gvar.sibry-gvar.simag);
%% 
[ne1,rho1]=calc_Ne(nedr,eqdsk_r,rhopsi,ne,rho);
[Ne]=calc_nete1(ne,ne1,rho,rhopsi,rho1);
figure(1)
hold on;
%contourf(eqdsk_r,eqdsk_z,Ne,50);
%contourf(eqdsk_r,eqdsk_z,-gvar.psirz',50);
%Ne=calc_nete(ne,rho,rhopsi);
Te=calc_nete(Te,rho,rhopsi);
%% 计算
for i=1:129
    for j=1:129
%         for k=1:101
            if rhopsi(i)>0.8 && rhopsi(i)<1.8
                
                deltan(i)=0.3;

            else
                deltan(i)=0;
            end  
%         end  
    end
end
%%


[Bt1]=calc_interp2(Bt);
[Br1]=calc_interp2(Br);
[Bz1]=calc_interp2(Bz);
[Ne1]=calc_interp2(Ne);
[Te1]=calc_interp2(Te);
%[deltan1]=calc_interp2(deltan);
[dpsidr]=calc_interp2(dpsidr);
[dpsidz]=calc_interp2(dpsidz);


Btot1=(Br1.^2+Bz1.^2+Bt1.^2).^0.5;
Ni1=Ne1/1.25;
[rhopsi11]=calc_interp2(rhopsi);
% [Ne]=calc_interp(Ne);



%global Bt1 Br1 Bz1 eqdsk_r1 eqdsk_z1  Ne1 dNe1 dBrr1 dBzz1 w Btot1
global Bt1 Br1 Bz1 eqdsk_r1 eqdsk_z1  Ne1 Te1 w Btot1 Ni1 deltan1
w=w1;

%[Ne]=calc_Ne(nedr);
%dNe1=gradient(Ne1);
%%

parfor j=1:nray

    npar=-2.04;
    a=0.45;
    ksi0roui=0.1;
    cirho=1;

    z=0;
    r=2.33;
    phi=0;
    
    [ii,jj]=findposition(z,r);
    Btot=[];

    nperps=[];
    nperpf=[];
    cnteta1=[];
    z1=[];
    r1=[];
    phi1=[];
    cnz1=[];
    cnr1=[];
    cm1=[];
    us1=[];
    Npar1=[];
    teta1=[];
    rhopsi1=[];
    
    dddd1=[];
    ddt1=[];
    ii1=[];
    jj1=[];
    sne1=[];
    npar_acc1=[];
    sbr=[];
    sbz=[];
    sbt=[];
    nperps1=[];
    nperpf1=[];
    nperp=[];
    
%%

%Btot(ii,jj)=(Br1(ii,jj)^2+Bz1(ii,jj)^2+Bt1(ii,jj)^2)^0.5;
[xar1,xar2,yar1,yar2,npar_acc]=calc_xyar(Ne1(ii,jj),Ni1(ii,jj),Btot1(ii,jj));
[eps]=calc_eps1(xar1,xar2,yar1,yar2);
[nperps,nperpf,D]=calc_nperpsf(eps,npar);
%[eps,nperps(1),nperpf(1),D]=calc_eps(xar1(1),xar2(1),yar1(1),yar2(1),npar);
kpar=w*npar./c;

if imag(nperps)==0
    nperp=nperps;
end
 gradpsi=sqrt(dpsidr.^2+dpsidz.^2);
% B_theta=(Bz.*dpsidr-Br.*dpsidz)./gradpsi;
cnteta1(1)=0;
cnphi=npar*Btot1(ii,jj)/Bt1(ii,jj);
cn2=npar^2+nperp^2;
cnrho2=cn2-cnteta1(1)^2-cnphi^2;
cnrho=sqrt(cnrho2);
cnr=(cirho*cnrho*dpsidr(ii,jj))/gradpsi(ii,jj);
cnz=(cirho*cnrho*dpsidz(ii,jj))/gradpsi(ii,jj);
cm=cnphi*eqdsk_r1(jj);

npar1=(cnr*Br1(ii,jj)+cnz*Bz1(ii,jj)+Bt1(ii,jj)*cnphi)/Btot1(ii,jj);
us=0;
z1(1)=z;
r1(1)=r;
phi1(1)=phi;
cnz1(1)=cnz;
cnr1(1)=cnr;
cm1(1)=cm;
us1(1)=us;
Npar1(1)=npar1;
teta1(1)=0;


%%
for i=1:lengthray
%     complex(j,i)
   
%     [ddd,ddt,ii,jj]=ddd1(z,r,phi,cnz,cnr,cm);
%     dddd(i)=ddd;
%     ddt1(i)=ddt;
%     ii1(i)=ii;
%     jj1(i)=jj;
%     sne1(i)=Ne1(ii,jj);
%     Btot(i)=sqrt(Br1(ii,jj)^2+Bz1(ii,jj)^2+Bt1(ii,jj)^2);
    
%%
%if (Ne1(ii,jj)==0 || r<eqdsk_r(1) || r>eqdsk_r(end) || ddt.z>c || ddt.r>c) %ddt>c 属于什么问题？
if inpolygon(r,z,gvar.rlim,gvar.zlim) && inpolygon(r,z,gvar.rbbbs,gvar.zbbbs)==0
%if inpolygon(r,z,gvar.rlim,gvar.zlim)
    [ddd,ddt,ii,jj]=calc_ddd1(z,r,phi,cnz,cnr,cm);
    if abs(ddt.z)<c && abs(ddt.r)<c
%     dddd1(i)=ddd;
%     ddt1(i)=ddt;
    ii1(i)=ii;
    jj1(i)=jj;
    sne1(i)=Ne1(ii,jj);
    Btot(i)=sqrt(Br1(ii,jj)^2+Bz1(ii,jj)^2+Bt1(ii,jj)^2);
    [xar1,xar2,yar1,yar2,npar_acc]=calc_xyar(sne1(i),sne1(i),Btot(i));
    [eps]=calc_eps1(xar1,xar2,yar1,yar2);
    [nperps,nperpf,D]=calc_nperpsf(eps,Npar1(i));
    npar_acc1(i)=npar_acc;
    sbr(i)=Br1(ii,jj);
    sbz(i)=Bz1(ii,jj);
    sbt(i)=Bt1(ii,jj);
    nperps1(i)=imag(nperps);
    nperpf1(i)=imag(nperpf);
%%
[z,r,phi,cnz,cnr,cm,us]=rk2(z,r,phi,cnz,cnr,cm,us);
[ii,jj]=findposition(z,r);
z1(i+1)=z;
r1(i+1)=r;
phi1(i+1)=phi;
cnz1(i+1)=cnz;
cnr1(i+1)=cnr;
cm1(i+1)=cm;
us1(i+1)=us;
[teta]=calc_teta11(z1(i+1),r1(i+1),z1(i),r1(i),gvar);
teta1(i+1)=teta1(i)+teta;
% if teta1(i)>teta
%     teta1(i+1)=teta+2*pi;
% else
%     teta1(i+1)=teta;
% end

Npar1(i+1)=(Br1(ii,jj)*cnr+Bz1(ii,jj)*cnz+Bt1(ii,jj)*cm/eqdsk_r1(jj))/Btot(i);
% z1(i+1)=z+tstep*ddt.z;
% r1(i+1)=r+tstep*ddt.r;
% phi1(i+1)=phi+tstep*ddt.phi;
% cnz1(i+1)=cnz+tstep*ddt.cnz;
% cnr1(i+1)=cnr+tstep*ddt.cnr;
% cm1(i+1)=cm+tstep*ddt.cm;
% us1(i+1)=us+sqrt((tstep*ddt.z)^2+(tstep*ddt.r)^2+(tstep*ddt.phi)^2);
% z=z1(i+1);
% r=r1(i+1);
% phi=phi1(i+1);
% cnz=cnz1(i+1);
% cnr=cnr1(i+1);
% cm=cm1(i+1);
% us=us1(i+1);

%%  区分刮削层和内部涨落水平
% if  j~=1
%     deltan11(i)=deltan1(ii,jj);
%     if  deltan1(ii,jj)~=0  &&  (imag(nperps)~=0 && imag(nperpf)~=0)~=1
%         [snperm,ds,cnz,cnr,cm,xrandom]=densfluc(z,r,phi,cnz,cnr,cm,tstep);       %在这里并没有给出新的（cnz,cnr,cm）,已修改
%         xrandom1(i)=xrandom;
%         snperm1(i)=snperm;
%        
%     end
% else
%     ds=0;
%end
%  us=us+ds;      %意义在哪？如果(z,r)不变的话，那么ds的变化是怎么来的？

%% 
%[z,r,phi,cnz,cnr,cm,us]=rk11(z,r,phi,cnz,cnr,cm,us);
%[z,r,phi,cnz,cnr,cm,us]=rk2(z,r,phi,cnz,cnr,cm,us);



%     if ((abs(Npar1(i))<npar_acc1(i))||(imag(nperps)>0&&imag(nperpf)>0))     %截止条件
% %         Npar1=Npar(i);
% %         z1(i)=z;
% %         r1(i)=r;
% %         phi1(i)=phi;
% %         us1(i)=us;
% %         cnz1(i)=cnz;
% %         cnr1(i)=cnr;
% %         cm1(i)=cm;
% %         rhopsi1(i)=rhopsi11(ii,jj);
% %         cnteta(i)=(cnr*Br1(ii,jj)+cnz*Bz1(ii,jj))/sqrt(Br1(ii,jj)^2+Bz1(ii,jj)^2);
%         rhopsi1(i+1)=rhopsi11(ii,jj);
%         cnteta1(i+1)=(cnr*Br1(ii,jj)+cnz*Bz1(ii,jj))/sqrt(Br1(ii,jj)^2+Bz1(ii,jj)^2);
%         break;
%     else
%         if rhopsi11(ii,jj)>1.2      %在rho为1.2处发生反射
%             cnteta=(cnr*Br1(ii,jj)+cnz*Bz1(ii,jj))/sqrt(Br1(ii,jj)^2+Bz1(ii,jj)^2);
%             cnpsi=(cnz*Br1(ii,jj)-cnr*Bz1(ii,jj))/sqrt(Br1(ii,jj)^2+Bz1(ii,jj)^2);
%             cnzref=(cnteta*Bz1(ii,jj)-cnpsi*Br1(ii,jj))/sqrt(Br1(ii,jj)^2+Bz1(ii,jj)^2);
%             cnrref=(cnteta*Br1(ii,jj)+cnpsi*Bz1(ii,jj))/sqrt(Br1(ii,jj)^2+Bz1(ii,jj)^2);
%             cnz=cnzref;
%             cnr=cnrref;
% %         end
         rhopsi1(i+1)=rhopsi11(ii,jj);
         cnteta1(i+1)=(cnr*Br1(ii,jj)+cnz*Bz1(ii,jj))/sqrt(Br1(ii,jj)^2+Bz1(ii,jj)^2);
%     end
    end
end
%%
end
end
%plot1(r1,z1,us1,Npar,cnteta,j);