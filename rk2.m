function [z,r,phi,cnz,cnr,cm,us]=rk2(z,r,phi,cnz,cnr,cm,us)
% t=rk(1);
% tend=rk(2);
% step=rk(3);
% eps=rk(4);
step=1e-9;


alpha= [1.0/4.0 3.0/8.0 12.0/13.0 1.0 1.0/2.0];
beta=[1.0/4.0, 3.0/32.0,  1932.0/2197.0,   8341.0/4104.0,  -6080.0/20520.0;
      0.0, 9.0/32.0, -7200.0/2197.0, -32832.0/4104.0,  41040.0/20520.0;
      0.0,      0.0,  7296.0/2197.0,  29440.0/4104.0, -28352.0/20520.0;
      0.0,      0.0,            0.0,   -845.0/4104.0,   9295.0/20520.0;
      0.0,      0.0,            0.0,             0.0,  -5643.0/20520.0 ];

  gama =[  902880.0/7618050.0,   -2090.0/752400.0;
                      0.0,                0.0;
      3953664.0/7618050.0,   22528.0/752400.0;
      3855735.0/7618050.0,   21970.0/752400.0;
     -1371249.0/7618050.0,  -15048.0/752400.0;
       277020.0/7618050.0,  -27360.0/752400.0];

zt=z;
rt=r;
phit=phi;
cnzt=cnz;
cnrt=cnr;
cmt=cm;
[~,ddt1,ii,jj]=calc_ddd1(zt,rt,phit,cnzt,cnrt,cmt);
zt=z+step*ddt1.z*beta(1,1);
rt=r+step*ddt1.r*beta(1,1);
phit=phi+step*ddt1.phi*beta(1,1);
cnzt=cnz+step*ddt1.cnz*beta(1,1);
cnrt=cnr+step*ddt1.cnr*beta(1,1);
cmt=cm+step*ddt1.cm*beta(1,1);
[~,ddt2,ii,jj]=calc_ddd1(zt,rt,phit,cnzt,cnrt,cmt);
ddt2.z=ddt2.z*step;
ddt2.r=ddt2.r*step;
ddt2.phi=ddt2.phi*step;
ddt2.cnz=ddt2.cnz*step;
ddt2.cnr=ddt2.cnr*step;
ddt2.cm=ddt2.cm*step;

zt=z+step*(ddt1.z*beta(1,2)+ddt2.z*beta(2,2));
rt=r+step*(ddt1.r*beta(1,2)+ddt2.r*beta(2,2));
phit=phi+step*(ddt1.phi*beta(1,2)+ddt2.phi*beta(2,2));
cnzt=cnz+step*(ddt1.cnz*beta(1,2)+ddt2.cnz*beta(2,2));
cnrt=cnr+step*(ddt1.cnr*beta(1,2)+ddt2.cnr*beta(2,2));
cmt=cm+step*(ddt1.cm*beta(1,2)+ddt2.cm*beta(2,2));
[~,ddt3,ii,jj]=calc_ddd1(zt,rt,phit,cnzt,cnrt,cmt);
%ddt3=ddt3*step;
ddt3.z=ddt3.z*step;
ddt3.r=ddt3.r*step;
ddt3.phi=ddt3.phi*step;
ddt3.cnz=ddt3.cnz*step;
ddt3.cnr=ddt3.cnr*step;
ddt3.cm=ddt3.cm*step;
zt=z+step*(ddt1.z*beta(1,3)+ddt2.z*beta(2,3)+ddt3.z*beta(3,3));
rt=r+step*(ddt1.r*beta(1,3)+ddt2.r*beta(2,3)+ddt3.r*beta(3,3));
phit=phi+step*(ddt1.phi*beta(1,3)+ddt2.phi*beta(2,3)+ddt3.phi*beta(3,3));
cnzt=cnz+step*(ddt1.cnz*beta(1,3)+ddt2.cnz*beta(2,3)+ddt3.cnz*beta(3,3));
cnrt=cnr+step*(ddt1.cnr*beta(1,3)+ddt2.cnr*beta(2,3)+ddt3.cnr*beta(3,3));
cmt=cm+step*(ddt1.cm*beta(1,3)+ddt2.cm*beta(2,3)+ddt3.cm*beta(3,3));
[~,ddt4,~,~]=calc_ddd1(zt,rt,phit,cnzt,cnrt,cmt);
%ddt4=ddt4*step;
ddt4.z=ddt4.z*step;
ddt4.r=ddt4.r*step;
ddt4.phi=ddt4.phi*step;
ddt4.cnz=ddt4.cnz*step;
ddt4.cnr=ddt4.cnr*step;
ddt4.cm=ddt4.cm*step;
zt=z+step*(ddt1.z*beta(1,4)+ddt2.z*beta(2,4)+ddt3.z*beta(3,4)+ddt4.z*beta(4,4))
rt=r+step*(ddt1.r*beta(1,4)+ddt2.r*beta(2,4)+ddt3.r*beta(3,4)+ddt4.r*beta(4,4))
phit=phi+step*(ddt1.phi*beta(1,4)+ddt2.phi*beta(2,4)+ddt3.phi*beta(3,4)+ddt4.phi*beta(4,4))
cnzt=cnz+step*(ddt1.cnz*beta(1,4)+ddt2.cnz*beta(2,4)+ddt3.cnz*beta(3,4)+ddt4.cnz*beta(4,4))
cnrt=cnr+step*(ddt1.cnr*beta(1,4)+ddt2.cnr*beta(2,4)+ddt3.cnr*beta(3,4)+ddt4.cnr*beta(4,4))
cmt=cm+step*(ddt1.cm*beta(1,4)+ddt2.cm*beta(2,4)+ddt3.cm*beta(3,4)+ddt4.cm*beta(4,4))
[~,ddt5,~,~]=calc_ddd1(zt,rt,phit,cnzt,cnrt,cmt)
%ddt5=ddt5*step;
ddt5.z=ddt5.z*step;
ddt5.r=ddt5.r*step;
ddt5.phi=ddt5.phi*step;
ddt5.cnz=ddt5.cnz*step;
ddt5.cnr=ddt5.cnr*step;
ddt5.cm=ddt5.cm*step;
zt=z+step*(ddt1.z*beta(1,5)+ddt2.z*beta(2,5)+ddt3.z*beta(3,5)+ddt4.z*beta(4,5)+ddt5.z*beta(5,5));
rt=r+step*(ddt1.r*beta(1,5)+ddt2.r*beta(2,5)+ddt3.r*beta(3,5)+ddt4.r*beta(4,5)+ddt5.r*beta(5,5));
phit=phi+step*(ddt1.phi*beta(1,5)+ddt2.phi*beta(2,5)+ddt3.phi*beta(3,5)+ddt4.phi*beta(4,5)+ddt5.phi*beta(5,5));
cnzt=cnz+step*(ddt1.cnz*beta(1,5)+ddt2.cnz*beta(2,5)+ddt3.cnz*beta(3,5)+ddt4.cnz*beta(4,5)+ddt5.cnz*beta(5,5));
cnrt=cnr+step*(ddt1.cnr*beta(1,5)+ddt2.cnr*beta(2,5)+ddt3.cnr*beta(3,5)+ddt4.cnr*beta(4,5)+ddt5.cnr*beta(5,5));
cmt=cm+step*(ddt1.cm*beta(1,5)+ddt2.cm*beta(2,5)+ddt3.cm*beta(3,5)+ddt4.cm*beta(4,5)+ddt5.cm*beta(5,5));
[~,ddt6,~,~]=calc_ddd1(zt,rt,phit,cnzt,cnrt,cmt);
%ddt6=ddt6*step;
ddt6.z=ddt6.z*step;
ddt6.r=ddt6.r*step;
ddt6.phi=ddt6.phi*step;
ddt6.cnz=ddt6.cnz*step;
ddt6.cnr=ddt6.cnr*step;
ddt6.cm=ddt6.cm*step;
ydelz=step*(ddt1.z*gama(1,2)+ddt2.z*gama(2,2)+ddt3.z*gama(3,2)+ddt4.z*gama(4,2)+ddt5.z*gama(5,2)+ddt6.z*gama(6,2));


dzdt=ddt1.z*gama(1,1)+ddt2.z*gama(2,1)+ddt3.z*gama(3,1)+ddt4.z*gama(4,1)+ddt5.z*gama(5,1)+ddt6.z*gama(6,1);
drdt=ddt1.r*gama(1,1)+ddt2.r*gama(2,1)+ddt3.r*gama(3,1)+ddt4.r*gama(4,1)+ddt5.r*gama(5,1)+ddt6.r*gama(6,1);
dphidt=ddt1.phi*gama(1,1)+ddt2.phi*gama(2,1)+ddt3.phi*gama(3,1)+ddt4.phi*gama(4,1)+ddt5.phi*gama(5,1)+ddt6.phi*gama(6,1);
dcnzdt=ddt1.cnz*gama(1,1)+ddt2.cnz*gama(2,1)+ddt3.cnz*gama(3,1)+ddt4.cnz*gama(4,1)+ddt5.cnz*gama(5,1)+ddt6.cnz*gama(6,1);
dcnrdt=ddt1.cnr*gama(1,1)+ddt2.cnr*gama(2,1)+ddt3.cnr*gama(3,1)+ddt4.cnr*gama(4,1)+ddt5.cnr*gama(5,1)+ddt6.cnr*gama(6,1);
dcmdt=ddt1.cm*gama(1,1)+ddt2.cm*gama(2,1)+ddt3.cm*gama(3,1)+ddt4.cm*gama(4,1)+ddt5.cm*gama(5,1)+ddt6.cm*gama(6,1);
dwdt=0;
%Sdot = Sdot1*gama(1,1) + Sdot2*gama(2,1) + Sdot3*gama(3,1) + Sdot4*gama(4,1) + Sdot5*gama(5,1) + Sdot6*gama(6,1);
% dzdt = s1[0]*gama[0][0] + s2[0]*gama[1][0] + s3[0]*gama[2][0] + s4[0]*gama[3][0] + s5[0]*gama[4][0] + s6[0]*gama[5][0];
% thetadot = s1[1]*gama[0][0] + s2[1]*gama[1][0] + s3[1]*gama[2][0] + s4[1]*gama[3][0] + s5[1]*gama[4][0] + s6[1]*gama[5][0];
% zdot = s1[2]*gama[0][0] + s2[2]*gama[1][0] + s3[2]*gama[2][0] + s4[2]*gama[3][0] + s5[2]*gama[4][0] + s6[2]*gama[5][0];
% krhodot = s1[3]*gama[0][0] + s2[3]*gama[1][0] + s3[3]*gama[2][0] + s4[3]*gama[3][0] + s5[3]*gama[4][0] + s6[3]*gama[5][0];
% mdot = s1[4]*gama[0][0] + s2[4]*gama[1][0] + s3[4]*gama[2][0] + s4[4]*gama[3][0] + s5[4]*gama[4][0] + s6[4]*gama[5][0];
% kzdot = s1[5]*gama[0][0] + s2[5]*gama[1][0] + s3[5]*gama[2][0] + s4[5]*gama[3][0] + s5[5]*gama[4][0] + s6[5]*gama[5][0];
% omegadot = 0.0;

z=z+step*dzdt;
r=r+step*drdt;
phi=phi+step*dphidt;
cnz=cnz+step*dcnzdt;
cnr=cnr+step*dcnrdt;
cm=cm+step*dcmdt;
us=us+sqrt((step*dzdt)^2+(step*drdt)^2+(step*dphidt)^2);
% z1(i)=z+step*dzdt;
% r1(i)=r+step*drdt;
% phi1(i)=phi+step*dphidt;
% cnz1(i)=cnz+step*dcnzdt;
% cnr1(i)=cnr+step*dcnrdt;
% cm1(i)=cm+step*dcmdt;
% us1(i)=us+sqrt((step*dzdt)^2+(step*drdt)^2+(step*dphidt)^2);
% z=z1(i);
% r=r1(i);
% phi=phi1(i);
% cnz=cnz1(i);
% cnr=cnr1(i);
% cm=cm1(i);
% us=us1(i);
  