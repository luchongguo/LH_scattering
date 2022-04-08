% solve the dD/dcnr problem, in order to solve the ray tracing problem
function [ddd,ddt,ii,jj]=ddd1(z,r,phi,cnz,cnr,cm)
global Bt1 Br1 Bz1 eqdsk_r1 eqdsk_z1 Ni1 Ne1 dNe1 dBrr1 dBzz1 w
step=0.01;
c=299792458;

% z=ray.z;
% r=ray.r;
% phi=ray.phi;
% cnz=ray.cnz;
% cnr=ray.cnr;
% cm=ray.cm;
[ii,jj]=findposition(z,r);
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
ni=Ne1(ii,jj);
% dne=dNe1(ii,jj);
% dbrr=dBrr1(ii,jj);
% dbzz=dBzz1(ii,jj);

%dBtt=0;
btot=sqrt(bphi^2+bz^2+br^2);
[xar1,xar2,yar1,yar2]=calc_xyar(ne,ni,btot);
% eps.perp=1+xar1/(yar1)^2-xar2;
% eps.par=1-xar1-xar2;
% eps.xy=xar1/yar1;
%[eps,~]=calc_eps1(xar1,xar2,yar1,yar2);

%%          ddd.cnz
cnzp=cnz+step;
cnp=sqrt(cnzp^2+cnr^2+cm^2/r^2);
gg=cnzp*bz+cnr*br+bphi*cm/r;
btot=sqrt(bphi^2+bz^2+br^2);
arg=gg/(btot*cnp);
gammap=acos(arg);
%[D0p]=hamit(eps,gammap,cnp,ii1,jj1);
%[D0p]=calc_hamit1(gammap,cnp,ii1,jj1);
[D0p,~]=calc_hamit1(xar1,xar2,yar1,yar2,gammap,cnp,ii1,jj1);

cnzm=cnz-step;
cnm=sqrt(cnzm^2+cnr^2+cm^2/r^2);
gg=cnzm*bz+cnr*br+bphi*cm/r;
btot=sqrt(bphi^2+bz^2+br^2);
arg=gg/(btot*cnm);
gammam=acos(arg);
%[D0m]=hamit(eps,gammam,cnm,ii1,jj1);
[D0m,~]=calc_hamit1(xar1,xar2,yar1,yar2,gammam,cnm,ii1,jj1);
ddd.cnz=(D0p-D0m)/2/step;
%%          ddd.cnr
cnrp=cnr+step;
cnp=sqrt(cnrp^2+cnz^2+cm^2/r^2);
gg=cnz*bz+cnrp*br+bphi*cm/r;
btot=sqrt(bphi^2+bz^2+br^2);
arg=gg/(btot*cnp);
gammap=acos(arg);
%[D0p]=hamit(eps,gammap,cnp,ii1,jj1);
%[D0p]=calc_hamit1(gammap,cnp,ii1,jj1);
[D0p,~]=calc_hamit1(xar1,xar2,yar1,yar2,gammap,cnp,ii1,jj1);

cnrm=cnr-step;
cnm=sqrt(cnrm^2+cnz^2+cm^2/r^2);
gg=cnz*bz+cnrm*br+bphi*cm/r;
btot=sqrt(bphi^2+bz^2+br^2);
arg=gg/(btot*cnm);
gammam=acos(arg);
%[D0m]=hamit(eps,gammam,cnm,ii1,jj1);
[D0m,~]=calc_hamit1(xar1,xar2,yar1,yar2,gammam,cnm,ii1,jj1);
ddd.cnr=(D0p-D0m)/2/step;
%%          ddd.cm
cmp=cm+step;
cnp=sqrt(cnz^2+cnr^2+cmp^2/r^2);
gg=cnz*bz+cnr*br+bphi*cmp/r;
btot=sqrt(bphi^2+bz^2+br^2);
arg=gg/(btot*cnp);
gammap=acos(arg);
%[D0p]=hamit(eps,gammap,cnp,ii1,jj1);
%[D0p]=calc_hamit1(gammap,cnp,ii1,jj1);
[D0p,~]=calc_hamit1(xar1,xar2,yar1,yar2,gammap,cnp,ii1,jj1);

cmm=cm-step;
cnm=sqrt(cnz^2+cnr^2+cmm^2/r^2);
gg=cnz*bz+cnr*br+bphi*cmm/r;
btot=sqrt(bphi^2+bz^2+br^2);
arg=gg/(btot*cnm);
gammam=acos(arg);
%[D0m]=hamit(eps,gammam,cnm,ii1,jj1);
[D0m,~]=calc_hamit1(xar1,xar2,yar1,yar2,gammam,cnm,ii1,jj1);
ddd.cm=(D0p-D0m)/2/step;
%%      ddd.w
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
% eps.perp=1+xar1/(yar1)^2-xar2;
% eps.par=1-xar1-xar2;
% eps.xy=xar1/yar1;
%[eps,~]=calc_eps1(xar1,xar2,yar1,yar2);

cnp=sqrt(cnzp^2+cnrp^2+cmp^2/r^2);
gg=cnzp*bz+cnrp*br+bphi*cmp/r;
btot=sqrt(bphi^2+bz^2+br^2);
arg=gg/(btot*cnp);
gammap=acos(arg);
%[D0p]=hamit(eps,gammap,cnp,ii1,jj1);
[D0p,~]=calc_hamit1(xar1,xar2,yar1,yar2,gammap,cnp,ii1,jj1);

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
% eps.perp=1+xar1/(yar1)^2-xar2;
% eps.par=1-xar1-xar2;
% eps.xy=xar1/yar1;
%[eps,~]=calc_eps1(xar1,xar2,yar1,yar2);

cnm=sqrt(cnzm^2+cnrm^2+cmm^2/r^2);
gg=cnzm*bz+cnrm*br+bphi*cmm/r;
btot=sqrt(bphi^2+bz^2+br^2);
arg=gg/(btot*cnm);
gammam=acos(arg);
%[D0m]=hamit(eps,gammam,cnm,ii1,jj1);
[D0m,~]=calc_hamit1(xar1,xar2,yar1,yar2,gammam,cnm,ii1,jj1);
ddd.w=(D0p-D0m)/2/step/w;
% %%          ddd.r
% nep=ne+dne*step;
% nip=ni+dne*step;
% brp=br+dbrr*step;
% btotp=sqrt(bphi^2+bz^2+brp^2);
% [xar1,xar2,yar1,yar2]=calc_xyar(nep,nip,btotp);    
% % eps.perp=1+xar1/(yar1)^2-xar2;
% % eps.par=1-xar1-xar2;
% % eps.xy=xar1/yar1;
% [eps,~]=calc_eps1(xar1,xar2,yar1,yar2);
% 
% cnp=sqrt(cnz^2+cnr^2+cm^2/(rp^2);
% ggp=cnz*bz+cnrp*brp+bphi*cm/rp;
% 
% argp=ggp/(btotp*cnp);
% gammap=acos(argp);
% [D0p]=hamit(eps,gammap,cnp,ii1,jj1);
% 
% nem=ne-dne*step;
% nim=ni-dne*step;
% brm=br-dbrr*step;
% btotm=sqrt(bphi^2+bz^2+brm^2);
% [xar1,xar2,yar1,yar2]=calc_xyar(nem,nim,btotm);    %maybe btot is the lastest btot
% % eps.perp=1+xar1/(yar1)^2-xar2;
% % eps.par=1-xar1-xar2;
% % eps.xy=xar1/yar1;
% [eps,~]=calc_eps1(xar1,xar2,yar1,yar2);
% 
% cnm=sqrt(cnz^2+cnrm^2+cm^2/r^2);
% ggm=cnz*bz+cnrp*brm+bphi*cm/r;
% 
% argm=gg/(btotm*cnm);
% gammam=acos(argm);
% [D0m]=hamit(eps,gammam,cnm);
% ddd.r=(D0p-D0m)/2/step;

%%          ddd.r
rp=r+step;
[ii,jj]=findposition(z,rp);
brp=Br1(ii,jj);
bzp=Bz1(ii,jj);
bphip=Bt1(ii,jj);
nep=Ne1(ii,jj);
nip=Ne1(ii,jj);
btotp=sqrt(bphip^2+bzp^2+brp^2);
[xar1,xar2,yar1,yar2]=calc_xyar(nep,nip,btotp);    
%[eps,~]=calc_eps1(xar1,xar2,yar1,yar2);
cnp=sqrt(cnz^2+cnr^2+cm^2/(rp^2));
ggp=cnz*bzp+cnr*brp+bphip*cm/rp;
argp=ggp/(btotp*cnp);
gammap=acos(argp);
%[D0p]=hamit(eps,gammap,cnp,ii1,jj1);
[D0p,~]=calc_hamit1(xar1,xar2,yar1,yar2,gammap,cnp,ii,jj);

rm=r+step;
[ii,jj]=findposition(z,rm);
brm=Br1(ii,jj);
bzm=Bz1(ii,jj);
bphim=Bt1(ii,jj);
nem=Ne1(ii,jj);
nim=Ne1(ii,jj);
btotm=sqrt(bphim^2+bzm^2+brm^2);
[xar1,xar2,yar1,yar2]=calc_xyar(nem,nim,btotm);    %maybe btot is the lastest btot
%[eps,~]=calc_eps1(xar1,xar2,yar1,yar2);
cnm=sqrt(cnz^2+cnr^2+cm^2/rm^2);
ggm=cnz*bzm+cnr*brm+bphim*cm/r;

argm=ggm/(btotm*cnm);
gammam=acos(argm);
%[D0m]=hamit(eps,gammam,cnm,ii1,jj1);
[D0m,~]=calc_hamit1(xar1,xar2,yar1,yar2,gammam,cnm,ii,jj);
ddd.r=(D0p-D0m)/2/step;

%%          ddd.z
% nep=ne+dne*step;
% nip=ni+dne*step;
% bzp=bz+dbzz*step;
% btotp=sqrt(bphi^2+bzp^2+br^2);
% [xar1,xar2,yar1,yar2]=calc_xyar(nep,nip,btotp);    
% % eps.perp=1+xar1/(yar1)^2-xar2;
% % eps.par=1-xar1-xar2;
% % eps.xy=xar1/yar1;
% [eps,~]=calc_eps1(xar1,xar2,yar1,yar2);
% cnp=sqrt(cnzp^2+cnr^2+cm^2/r^2);
% ggp=cnzm*bzp+cnrp*br+bphi*cm/r;
% 
% argp=ggp/(btotp*cnp);
% gammap=acos(argp);
% [D0p]=hamit(eps,gammap,cnp,ii1,jj1);
% 
% nem=ne-dne*step;
% nim=ni-dne*step;
% bzm=bz-dbzz*step;
% btotm=sqrt(bphi^2+bzm^2+br^2);
% [xar1,xar2,yar1,yar2]=calc_xyar(nem,nim,btotm);    %maybe btot is the lastest btot
% % eps.perp=1+xar1/(yar1)^2-xar2;
% % eps.par=1-xar1-xar2;
% % eps.xy=xar1/yar1;
% [eps,~]=calc_eps1(xar1,xar2,yar1,yar2);
% cnm=sqrt(cnzm^2+cnr^2+cm^2/r^2);
% ggm=cnzm*bzm+cnrp*br+bphi*cm/r;
% 
% argm=ggm/(btotm*cnm);
% gammam=acos(argm);
% [D0m]=hamit(eps,gammam,cnm);
% ddd.z=(D0p-D0m)/2/step;
%%
zp=z+step;
[ii,jj]=findposition(zp,r);
brp=Br1(ii,jj);
bzp=Bz1(ii,jj);
bphip=Bt1(ii,jj);
nep=Ne1(ii,jj);
nip=Ne1(ii,jj);
btotp=sqrt(bphip^2+bzp^2+brp^2);
[xar1,xar2,yar1,yar2]=calc_xyar(nep,nip,btotp);    
%[eps,~]=calc_eps1(xar1,xar2,yar1,yar2);
cnp=sqrt(cnz^2+cnr^2+cm^2/(r^2));
ggp=cnz*bzp+cnr*brp+bphip*cm/r;
argp=ggp/(btotp*cnp);
gammap=acos(argp);
%[D0p]=hamit(eps,gammap,cnp,ii1,jj1);
[D0p,~]=calc_hamit1(xar1,xar2,yar1,yar2,gammap,cnp,ii,jj);

zm=z-step;
[ii,jj]=findposition(zm,r);
brm=Br1(ii,jj);
bzm=Bz1(ii,jj);
bphim=Bt1(ii,jj);
nem=Ne1(ii,jj);
nim=Ne1(ii,jj);
btotm=sqrt(bphim^2+bzm^2+brm^2);
[xar1,xar2,yar1,yar2]=calc_xyar(nem,nim,btotm);    %maybe btot is the lastest btot
%[eps,~]=calc_eps1(xar1,xar2,yar1,yar2);
cnm=sqrt(cnz^2+cnr^2+cm^2/r^2);
ggm=cnz*bzm+cnr*brm+bphim*cm/r;

argm=ggm/(btotm*cnm);
gammam=acos(argm);
%[D0m]=hamit(eps,gammam,cnm,ii1,jj1);
[D0m,~]=calc_hamit1(xar1,xar2,yar1,yar2,gammam,cnm,ii,jj);
ddd.z=(D0p-D0m)/2/step;
%%          ddd.phi
ddd.phi=0;
% ne=ne+dne*step;
% ni=ni+dne*step;
% bphi=bphi+dBtt*step;
% btot=sqrt(bphi^2+bz^2+br^2);
% [xar1,xar2,yar1,yar2]=calc_xyar(ne,ni,btot);    
% eps.perp=1+xar1/(yar1)^2-xar2;
% eps.par=1-xar1-xar2;
% eps.xy=xar1/yar1;
% cnp=sqrt(cnz^2+cnr^2+cm^2/r^2);
% gg=cnz*bz+cnrp*br+bphi*cm/r;
% arg=gg/(btot*cnp);
% gammap=acos(arg);
% [D0p]=hamit(eps,gammap,cnp,ii1,jj1);
% 
% ne=ne-dne*step;
% ni=ni-dne*step;
% bphi=bphi-dBtt*step;
% btot=sqrt(bphi^2+bz^2+br^2);
% [xar1,xar2,yar1,yar2]=calc_xyar(ne,ni,btot);    %maybe btot is the lastest btot
% eps.perp=1+xar1/(yar1)^2-xar2;
% eps.par=1-xar1-xar2;
% eps.xy=xar1/yar1;
% cnm=sqrt(cnz^2+cnr^2+cm^2/r^2);
% gg=cnz*bz+cnrp*br+bphi*cm/r;
% 
% arg=gg/(btot*cnm);
% gammam=acos(arg);
% [D0m]=hamit(eps,gammam,cnm);
% ddd.phi=(D0p-D0m)/2/step;
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


