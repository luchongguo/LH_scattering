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
nray=5;
lengthray=5000;
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
            if rhopsi(i,j)>0.8 && rhopsi(i,j)<1.8
                
                deltan(i,j)=0.3;

            else
                deltan(i,j)=0;
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
[deltan1]=calc_interp2(deltan);
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
%%

%Btot(ii,jj)=(Br1(ii,jj)^2+Bz1(ii,jj)^2+Bt1(ii,jj)^2)^0.5;
[xar1,xar2,yar1,yar2,npar_acc]=calc_xyar(Ne1(ii,jj),Ni1(ii,jj),Btot1(ii,jj));
[eps,reps]=calc_eps1(xar1,xar2,yar1,yar2);
[nperps,nperpf,D]=calc_nperpsf(eps,npar);
%[eps,nperps(1),nperpf(1),D]=calc_eps(xar1(1),xar2(1),yar1(1),yar2(1),npar);
kpar=w*npar./c;

if imag(nperps(1))==0
    nperp=nperps;
end
 gradpsi=sqrt(dpsidr.^2+dpsidz.^2);
% B_theta=(Bz.*dpsidr-Br.*dpsidz)./gradpsi;
cnteta1(1,j)=0;
cnphi=npar*Btot1(ii,jj)/Bt1(ii,jj);
cn2=npar^2+nperp^2;
cnrho2=cn2-cnteta1(1,j)^2-cnphi^2;
cnrho=sqrt(cnrho2);
cnr=(cirho*cnrho*dpsidr(ii,jj))/gradpsi(ii,jj);
cnz=(cirho*cnrho*dpsidz(ii,jj))/gradpsi(ii,jj);
cm=cnphi*eqdsk_r1(jj);

npar1=(cnr*Br1(ii,jj)+cnz*Bz1(ii,jj)+Bt1(ii,jj)*cnphi)/Btot1(ii,jj);
us=0;
z1(1,j)=z;
r1(1,j)=r;
phi1(1,j)=phi;
cnz1(1,j)=cnz;
cnr1(1,j)=cnr;
cm1(1,j)=cm;
us1(1,j)=us;
Npar1(1,j)=npar1;
teta1(1,j)=0;



%%
for i=1:lengthray
    complex(j,i)
   
%     [ddd,ddt,ii,jj]=ddd1(z,r,phi,cnz,cnr,cm);
%     dddd(i,j)=ddd;
%     ddt1(i,j)=ddt;
%     ii1(i,j)=ii;
%     jj1(i,j)=jj;
%     sne11(i,j)=Ne1(ii,jj);
%     Btot(i,j)=sqrt(Br1(ii,jj)^2+Bz1(ii,jj)^2+Bt1(ii,jj)^2);
    
%%
%if (Ne1(ii,jj)==0 || r<eqdsk_r(1) || r>eqdsk_r(end) || ddt.z>c || ddt.r>c) %ddt>c 属于什么问题？
if inpolygon(r,z,gvar.rlim,gvar.zlim) && inpolygon(r,z,gvar.rbbbs,gvar.zbbbs)==0
%if inpolygon(r,z,gvar.rlim,gvar.zlim)
    [ddd,ddt,ii,jj]=calc_ddd1(z,r,phi,cnz,cnr,cm);
    if abs(ddt.z)<c && abs(ddt.r)<c
    dddd1(i,j)=ddd;
    ddt1(i,j)=ddt;
    ii1(i,j)=ii;
    jj1(i,j)=jj;
    %sne11(i,j)=Ne1(ii,jj);
    Btot(i,j)=sqrt(Br1(ii,jj)^2+Bz1(ii,jj)^2+Bt1(ii,jj)^2);
    [xar1,xar2,yar1,yar2,npar_acc]=calc_xyar(Ne1(ii,jj),Ne1(ii,jj),Btot(i,j));
    [eps,reps]=calc_eps1(xar1,xar2,yar1,yar2);
    [nperps,nperpf,D]=calc_nperpsf(eps,tmpNpar1);
    %npar_acc1(i,j)=npar_acc;
    sbr(i,j)=Br1(ii,jj);
    sbz(i,j)=Bz1(ii,jj);
    sbt(i,j)=Bt1(ii,jj);
    nperps1(i,j)=imag(nperps);
    nperpf1(i,j)=imag(nperpf);
%%
[z1,r1,phi1,cnz1,cnr1,cm1,us1]=rk2(z,r,phi,cnz,cnr,cm,us);
[ii,jj]=findposition(z,r);
tmpz1=z1;
tmpr1=r1;
tmpphi1=phi1;
tmpcnz1=cnz1;
tmpcnr1=cnr1;
tmpcm1=cm1;
tmpus1=us1;
[teta1]=calc_teta11(tmpz1,tmpr1,z,r,gvar);
tmpteta1=tmpteta1+teta;
% if teta1(i,j)>teta
%     teta1(i+1,j)=teta+2*pi;
% else
%     teta1(i+1,j)=teta;
% end

tmpNpar1=(Br1(ii,jj)*cnr+Bz1(ii,jj)*cnz+Bt1(ii,jj)*cm/eqdsk_r1(jj))/Btot(i,j);
% tmpz1=z+tstep*ddt.z;
% tmpr1=r+tstep*ddt.r;
% tmpphi1=phi+tstep*ddt.phi;
% tmpcnz1=cnz+tstep*ddt.cnz;
% tmpcnr1=cnr+tstep*ddt.cnr;
% tmpcm1=cm+tstep*ddt.cm;
% tmpus1=us+sqrt((tstep*ddt.z)^2+(tstep*ddt.r)^2+(tstep*ddt.phi)^2);
% z=tmpz1;
% r=tmpr1;
% phi=tmpphi1;
% cnz=tmpcnz1;
% cnr=tmpcnr1;
% cm=tmpcm1;
% us=tmpus1;

%%  区分刮削层和内部涨落水平
if  j~=1
    deltan11(i,j)=deltan1(ii,jj);
    if  deltan1(ii,jj)~=0  &&  (imag(nperps)~=0 && imag(nperpf)~=0)~=1
        [snperm,ds,cnz,cnr,cm,xrandom]=densfluc(z,r,phi,cnz,cnr,cm,tstep);       %在这里并没有给出新的（cnz,cnr,cm）,已修改
        xrandom1(i,j)=xrandom;
        snperm1(i,j)=snperm;
       
    end
% else
%     ds=0;
end
%  us=us+ds;      %意义在哪？如果(z,r)不变的话，那么ds的变化是怎么来的？

%% 
%[z,r,phi,cnz,cnr,cm,us]=rk11(z,r,phi,cnz,cnr,cm,us);
%[z,r,phi,cnz,cnr,cm,us]=rk2(z,r,phi,cnz,cnr,cm,us);



    if ((abs(tmpNpar1)<npar_acc)||(imag(nperps)>0&&imag(nperpf)>0))     %截止条件
%         Npar1=Npar(i,j);
%         z1(i,j)=z;
%         r1(i,j)=r;
%         phi1(i,j)=phi;
%         us1(i,j)=us;
%         cnz1(i,j)=cnz;
%         cnr1(i,j)=cnr;
%         cm1(i,j)=cm;
%         rhopsi1(i,j)=rhopsi11(ii,jj);
%         cnteta(i,j)=(cnr*Br1(ii,jj)+cnz*Bz1(ii,jj))/sqrt(Br1(ii,jj)^2+Bz1(ii,jj)^2);
        tmprhopsi1=rhopsi11(ii,jj);
        tmpcnteta1=(cnr*Br1(ii,jj)+cnz*Bz1(ii,jj))/sqrt(Br1(ii,jj)^2+Bz1(ii,jj)^2);
        break;
    else
        if rhopsi11(ii,jj)>1.2      %在rho为1.2处发生反射
            cnteta=(cnr*Br1(ii,jj)+cnz*Bz1(ii,jj))/sqrt(Br1(ii,jj)^2+Bz1(ii,jj)^2);
            cnpsi=(cnz*Br1(ii,jj)-cnr*Bz1(ii,jj))/sqrt(Br1(ii,jj)^2+Bz1(ii,jj)^2);
            cnzref=(cnteta*Bz1(ii,jj)-cnpsi*Br1(ii,jj))/sqrt(Br1(ii,jj)^2+Bz1(ii,jj)^2);
            cnrref=(cnteta*Br1(ii,jj)+cnpsi*Bz1(ii,jj))/sqrt(Br1(ii,jj)^2+Bz1(ii,jj)^2);
            cnz=cnzref;
            cnr=cnrref;
        end
        tmprhopsi1=rhopsi11(ii,jj);
        tmpcnteta1=(cnr*Br1(ii,jj)+cnz*Bz1(ii,jj))/sqrt(Br1(ii,jj)^2+Bz1(ii,jj)^2);
    end
    end
end
%%
end
end
%plot1(r1,z1,us1,Npar,cnteta,j);