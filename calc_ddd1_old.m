% solve the dD/dcnr problem, in order to solve the ray tracing problem
function [ddd,ddt,ii,jj]=calc_ddd1_old(z,r,phi,cnz,cnr,cm)
%global Bt Bp Br Bz eqdsk_r eqdsk_z  Ne Ni Te w Btot Deltan rhopsi
global Bt1 Bp1 Br1 Bz1 eqdsk_r1 eqdsk_z1  Ne1 Ni1 Te1 w Btot1 Deltan1
step=0.001;  %如何确定这个量
c=299792458;

% z=rayz;
% r=rayr;
% phi=rayphi;
% cnz=raycnz;
% cnr=raycnr;
% cm=raycm;
%[ii,jj]=findposition(z,r);
[~,ii]=min(abs(z-eqdsk_z1));
[~,jj]=min(abs(r-eqdsk_r1));
% [imax,imin,jmax,jmin]=calc_position(z,r,eqdsk_z,eqdsk_r);
% calc_position_new;
ii1=ii;
jj1=jj;
% if ii==0||jj==0
%     ddd=0;
%     ddt=0;
% end

br=Br1(ii,jj);
bz=Bz1(ii,jj);
bphi=Bt1(ii,jj);
ne=Ne1(ii,jj);
ni=Ni1(ii,jj);
% dne=dNe1(ii,jj);
% dbrr=dBrr1(ii,jj);
% dbzz=dBzz1(ii,jj);

%dBtt=0;
btot=sqrt(bphi^2+bz^2+br^2);
[xar1,xar2,yar1,yar2]=calc_xyar(ne,ni,btot);
% epsperp=1+xar1/(yar1)^2-xar2;
% epspar=1-xar1-xar2;
% epsxy=xar1/yar1;
%[eps]=calc_eps1(xar1,xar2,yar1,yar2);

%%          dddcnz
cnzp=cnz+step;
[gammap,cnp]=calc_gamma(bz,br,bphi,cnzp,cnr,cm,r);
%[D0p]=hamit(eps,gammap,cnp,ii1,jj1);
%[D0p]=calc_D0_old(gammap,cnp,ii1,jj1);
[D0p]=calc_D0_old(xar1,xar2,yar1,yar2,gammap,cnp,ii1,jj1);

cnzm=cnz-step;
[gammam,cnm]=calc_gamma(bz,br,bphi,cnzm,cnr,cm,r);
%[D0m]=hamit(eps,gammam,cnm,ii1,jj1);
[D0m]=calc_D0_old(xar1,xar2,yar1,yar2,gammam,cnm,ii1,jj1);
ddd.cnz=(D0p-D0m)/2/step;
%%          dddcnr
cnrp=cnr+step;
[gammap,cnp]=calc_gamma(bz,br,bphi,cnz,cnrp,cm,r);
%[D0p]=hamit(eps,gammap,cnp,ii1,jj1);
%[D0p]=calc_D0_old(gammap,cnp,ii1,jj1);
[D0p]=calc_D0_old(xar1,xar2,yar1,yar2,gammap,cnp,ii1,jj1);

cnrm=cnr-step;
[gammam,cnm]=calc_gamma(bz,br,bphi,cnz,cnrm,cm,r);
%[D0m]=hamit(eps,gammam,cnm,ii1,jj1);
[D0m]=calc_D0_old(xar1,xar2,yar1,yar2,gammam,cnm,ii1,jj1);
ddd.cnr=(D0p-D0m)/2/step;
%%          dddcm
cmp=cm+step;
[gammap,cnp]=calc_gamma(bz,br,bphi,cnz,cnr,cmp,r);
%[D0p]=hamit(eps,gammap,cnp,ii1,jj1);
%[D0p]=calc_D0_old(gammap,cnp,ii1,jj1);
[D0p]=calc_D0_old(xar1,xar2,yar1,yar2,gammap,cnp,ii1,jj1);

cmm=cm-step;
[gammam,cnm]=calc_gamma(bz,br,bphi,cnz,cnr,cmm,r);
%[D0m]=hamit(eps,gammam,cnm,ii1,jj1);
[D0m]=calc_D0_old(xar1,xar2,yar1,yar2,gammam,cnm,ii1,jj1);
ddd.cm=(D0p-D0m)/2/step;
%%      dddw
wp=w+step*w;
dwp=w/wp;
[xar1,xar2,yar1,yar2]=calc_xyar(ne,ni,btot);
xar1=xar1*dwp*dwp;
xar2=xar2*dwp*dwp;
yar1=yar1*dwp;
yar2=yar2*dwp;
cnrp=cnr*dwp;
cnzp=cnz*dwp;
cmp=cm*dwp;
% epsperp=1+xar1/(yar1)^2-xar2;
% epspar=1-xar1-xar2;
% epsxy=xar1/yar1;
%[eps]=calc_eps1(xar1,xar2,yar1,yar2);

[gammap,cnp]=calc_gamma(bz,br,bphi,cnzp,cnrp,cmp,r);
%[D0p]=hamit(eps,gammap,cnp,ii1,jj1);
[D0p]=calc_D0_old(xar1,xar2,yar1,yar2,gammap,cnp,ii1,jj1);

wm=w-step*w;
dwm=w/wm;
[xar1,xar2,yar1,yar2]=calc_xyar(ne,ni,btot);
xar1=xar1*dwm*dwm;
xar2=xar2*dwm*dwm;
yar1=yar1*dwm;
yar2=yar2*dwm;
cnrm=cnr*dwm;
cnzm=cnz*dwm;
cmm=cm*dwm;
% epsperp=1+xar1/(yar1)^2-xar2;
% epspar=1-xar1-xar2;
% epsxy=xar1/yar1;
%[eps]=calc_eps1(xar1,xar2,yar1,yar2);

[gammam,cnm]=calc_gamma(bz,br,bphi,cnzm,cnrm,cmm,r);
%[D0m]=hamit(eps,gammam,cnm,ii1,jj1);
[D0m]=calc_D0_old(xar1,xar2,yar1,yar2,gammam,cnm,ii1,jj1);
ddd.w=(D0p-D0m)/2/step/w;

%%          dddr
rp=r+step;
%[ii,jj]=findposition(z,rp);
[~,ii]=min(abs(z-eqdsk_z1));
[~,jj]=min(abs(rp-eqdsk_r1));
brp=Br1(ii,jj);
bzp=Bz1(ii,jj);
bphip=Bt1(ii,jj);
nep=Ne1(ii,jj);
nip=Ni1(ii,jj);
btotp=sqrt(bphip^2+bzp^2+brp^2);
[xar1,xar2,yar1,yar2]=calc_xyar(nep,nip,btotp);    
%[eps]=calc_eps1(xar1,xar2,yar1,yar2);
[gammap,cnp]=calc_gamma(bzp,brp,bphip,cnz,cnr,cm,rp);
%[D0p]=hamit(eps,gammap,cnp,ii1,jj1);
[D0p]=calc_D0_old(xar1,xar2,yar1,yar2,gammap,cnp,ii,jj);

rm=r-step;
% [ii,jj]=findposition(z,rm);
[~,ii]=min(abs(z-eqdsk_z1));
[~,jj]=min(abs(rm-eqdsk_r1));
brm=Br1(ii,jj);
bzm=Bz1(ii,jj);
bphim=Bt1(ii,jj);
nem=Ne1(ii,jj);
nim=Ni1(ii,jj);
btotm=sqrt(bphim^2+bzm^2+brm^2);
[xar1,xar2,yar1,yar2]=calc_xyar(nem,nim,btotm);    %maybe btot is the lastest btot
%[eps]=calc_eps1(xar1,xar2,yar1,yar2);
[gammam,cnm]=calc_gamma(bzm,brm,bphim,cnz,cnr,cm,rm);
%[D0m]=hamit(eps,gammam,cnm,ii1,jj1);
[D0m]=calc_D0_old(xar1,xar2,yar1,yar2,gammam,cnm,ii,jj);
ddd.r=(D0p-D0m)/2/step;


%%
zp=z+step;
%[ii,jj]=findposition(zp,r);
[~,ii]=min(abs(zp-eqdsk_z1));
[~,jj]=min(abs(r-eqdsk_r1));
brp=Br1(ii,jj);
bzp=Bz1(ii,jj);
bphip=Bt1(ii,jj);
nep=Ne1(ii,jj);
nip=Ni1(ii,jj);
btotp=sqrt(bphip^2+bzp^2+brp^2);
[xar1,xar2,yar1,yar2]=calc_xyar(nep,nip,btotp);    
%[eps]=calc_eps1(xar1,xar2,yar1,yar2);
[gammap,cnp]=calc_gamma(bzp,brp,bphip,cnz,cnr,cm,rp);
%[D0p]=hamit(eps,gammap,cnp,ii1,jj1);
[D0p]=calc_D0_old(xar1,xar2,yar1,yar2,gammap,cnp,ii,jj);

zm=z-step;
%[ii,jj]=findposition(zm,r);
[~,ii]=min(abs(zm-eqdsk_z1));
[~,jj]=min(abs(r-eqdsk_r1));
brm=Br1(ii,jj);
bzm=Bz1(ii,jj);
bphim=Bt1(ii,jj);
nem=Ne1(ii,jj);
nim=Ni1(ii,jj);
btotm=sqrt(bphim^2+bzm^2+brm^2);
[xar1,xar2,yar1,yar2]=calc_xyar(nem,nim,btotm);    %maybe btot is the lastest btot
%[eps]=calc_eps1(xar1,xar2,yar1,yar2);
[gammam,cnm]=calc_gamma(bzm,brm,bphim,cnz,cnr,cm,r);
%[D0m]=hamit(eps,gammam,cnm,ii1,jj1);
[D0m]=calc_D0_old(xar1,xar2,yar1,yar2,gammam,cnm,ii,jj);
ddd.z=(D0p-D0m)/2/step;
%%          dddphi
phip=phi+step;
%[ii,jj]=findposition(z,r);
[~,ii]=min(abs(z-eqdsk_z1));
[~,jj]=min(abs(r-eqdsk_r1));
brp=Br1(ii,jj);
bzp=Bz1(ii,jj);
bphip=Bt1(ii,jj);
nep=Ne1(ii,jj);
nip=Ni1(ii,jj);
btotp=sqrt(bphip^2+bzp^2+brp^2);
[xar1,xar2,yar1,yar2]=calc_xyar(nep,nip,btotp);    
%[eps]=calc_eps1(xar1,xar2,yar1,yar2);
[gammap,cnp]=calc_gamma(bzp,brp,bphip,cnz,cnr,cm,r);
%[D0p]=hamit(eps,gammap,cnp,ii1,jj1);
[D0p]=calc_D0_old(xar1,xar2,yar1,yar2,gammap,cnp,ii,jj);

phi=phi-step;
%[ii,jj]=findposition(z,r);
[~,ii]=min(abs(z-eqdsk_z1));
[~,jj]=min(abs(r-eqdsk_r1));
brm=Br1(ii,jj);
bzm=Bz1(ii,jj);
bphim=Bt1(ii,jj);
nem=Ne1(ii,jj);
nim=Ni1(ii,jj);
btotm=sqrt(bphim^2+bzm^2+brm^2);
[xar1,xar2,yar1,yar2]=calc_xyar(nem,nim,btotm);    %maybe btot is the lastest btot
%[eps]=calc_eps1(xar1,xar2,yar1,yar2);
[gammam,cnm]=calc_gamma(bzm,brm,bphim,cnz,cnr,cm,r);
%[D0m]=hamit(eps,gammam,cnm,ii1,jj1);
[D0m]=calc_D0_old(xar1,xar2,yar1,yar2,gammam,cnm,ii,jj);
ddd.phi=(D0p-D0m)/2/step;

%%

%%
ddt.z=-ddd.cnz/ddd.w*c/w;
ddt.r=-ddd.cnr/ddd.w*c/w;
ddt.phi=-ddd.cm/ddd.w*c/w;
ddt.cnz=ddd.z/ddd.w*c/w;
ddt.cnr=ddd.r/ddd.w*c/w;
ddt.cm=ddd.phi/ddd.w*c/w;

ii=ii1;
jj=jj1;