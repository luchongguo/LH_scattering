function [z,r,phi,cnz,cnr,cm,us]=rk1(z,r,phi,cnz,cnr,cm,us)
% t=rk(1);
% tend=rk(2);
% step=rk(3);
% eps=rk(4);
step=1e-9; %不改变计算精度


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
[~,ddtz1]=ddd1(zt,r,phi,cnz,cnr,cm);
zt=z+step*ddtz1.z*beta(1,1);
[~,ddtz2]=ddd1(zt,r,phi,cnz,cnr,cm);
ddtz2=step*ddtz2;
zt=z+step*(ddtz1.z*beta(1,2)+ddtz2.z*beta(2,2));
[~,ddtz3]=ddd1(zt,r,phi,cnz,cnr,cm);
ddtz3=step*ddtz3;
zt=z+step*(ddtz1.z*beta(1,3)+ddtz2.z*beta(2,3)+ddtz3.z*beta(3,3));
[~,ddtz4]=ddd1(zt,r,phi,cnz,cnr,cm);
ddtz4=step*ddtz4;
zt=z+step*(ddtz1.z*beta(1,4)+ddtz2.z*beta(2,4)+ddtz3.z*beta(3,4)+ddtz4.z*beta(4,4));
[~,ddtz5]=ddd1(zt,r,phi,cnz,cnr,cm);
ddtz5=step*ddtz5;
zt=z+step*(ddtz1.z*beta(1,5)+ddtz2.z*beta(2,5)+ddtz3.z*beta(3,5)+ddtz4.z*beta(4,5)+ddtz5.z*beta(5,5));
[~,ddtz6]=ddd1(zt,r,phi,cnz,cnr,cm);
ddtz6=step*ddtz6;
ydelz=step*(ddtz1.z*gama(1,2)+ddtz2.z*gama(2,2)+ddtz3.z*gama(3,2)+ddtz4.z*gama(4,2)+ddtz5.z*gama(5,2)+ddtz6.z*gama(6,2));





% zt=z;
% [~,ddt1]=ddd1(zt,r,phi,cnz,cnr,cm);
% zt=z+step*ddt1.z*beta(1,1);
% [~,ddt2]=ddd1(zt,r,phi,cnz,cnr,cm);
% zt=z+step*(ddt1.z*beta(1,2)+ddt2.z*beta(2,2));
% [~,ddt3]=ddd1(zt,r,phi,cnz,cnr,cm);
% zt=z+step*(ddt1.z*beta(1,3)+ddt2.z*beta(2,3)+ddt3.z*beta(3,3));
% [~,ddt4]=ddd1(zt,r,phi,cnz,cnr,cm);
% zt=z+step*(ddt1.z*beta(1,4)+ddt2.z*beta(2,4)+ddt3.z*beta(3, 4)+ddt4.z*beta(4,4));
% [~,ddt5]=ddd1(zt,r,phi,cnz,cnr,cm);
% zt=z+step*(ddt1.z*beta(1,5)+ddt2.z*beta(2,5)+ddt3.z*beta(3,5)+ddt4.z*beta(4,5)+ddt5.z*beta(5,5));
% [~,ddt6]=ddd1(zt,r,phi,cnz,cnr,cm);
% ydelz=step*(ddt1.z*gama(1,2)+ddt2.z*gama(2,2)+ddt3.z*gama(3,2)+ddt4.z*gama(4,2)+ddt5.z*gama(5,2)+ddt6.z*gama(6,2));

% rt=r;
% [~,ddt1]=ddd1(z,rt,phi,cnz,cnr,cm);
% rt=z+step*ddt1.r*beta(1,1);
% [~,ddt2]=ddd1(z,rt,phi,cnz,cnr,cm);
% rt=z+step*(ddt1.r*beta(1,2)+ddt2.r*beta(2,2));
% [~,ddt3]=ddd1(z,rt,phi,cnz,cnr,cm);
% rt=z+step*(ddt1.r*beta(1,3)+ddt2.r*beta(2,3)+ddt3.r*beta(3,3));
% [~,ddt4]=ddd1(z,rt,phi,cnz,cnr,cm);
% rt=z+step*(ddt1.r*beta(1,4)+ddt2.r*beta(2,4)+ddt3.r*beta(3,4)+ddt4.r*beta(4,4));
% [~,ddt5]=ddd1(z,rt,phi,cnz,cnr,cm);
% rt=z+step*(ddt1.r*beta(1,5)+ddt2.r*beta(2,5)+ddt3.r*beta(3,5)+ddt4.r*beta(4,5)+ddt5.r*beta(5,5));
% [~,ddt6]=ddd1(z,rt,phi,cnz,cnr,cm);
% ydelr=step*(ddt1.r*gama(1,2)+ddt2.r*gama(2,2)+ddt3.r*gama(3,2)+ddt4.r*gama(4,2)+ddt5.r*gama(5,2)+ddt6.r*gama(6,2));

phit=phi;
[~,ddt1]=ddd1(z,r,phit,cnz,cnr,cm);
phit=phi+step*ddt1.phi*beta(1,1);
[~,ddt2]=ddd1(z,r,phit,cnz,cnr,cm);
ddt2=step*ddt2;
phit=phi+step*(ddt1.phi*beta(1,2)+ddt2.phi*beta(2,2));
[~,ddt3]=ddd1(z,r,phit,cnz,cnr,cm);
ddt3=step*ddt3;
phit=phi+step*(ddt1.phi*beta(1,3)+ddt2.phi*beta(2,3)+ddt3.phi*beta(3,3));
[~,ddt4]=ddd1(z,r,phit,cnz,cnr,cm);
ddt4=step*ddt4;
phit=phi+step*(ddt1.phi*beta(1,4)+ddt2.phi*beta(2,4)+ddt3.phi*beta(3,4)+ddt4.phi*beta(4,4));
[~,ddt5]=ddd1(z,r,phit,cnz,cnr,cm);
ddt5=step*ddt5;
phit=phi+step*(ddt1.phi*beta(1,5)+ddt2.phi*beta(2,5)+ddt3.phi*beta(3,5)+ddt4.phi*beta(4,5)+ddt5.phi*beta(5,5));
[~,ddt6]=ddd1(z,r,phit,cnz,cnr,cm);
ddt6=step*ddt6;
ydelphi=step*(ddt1.phi*gama(1,2)+ddt2.phi*gama(2,2)+ddt3.phi*gama(3,2)+ddt4.phi*gama(4,2)+ddt5.phi*gama(5,2)+ddt6.phi*gama(6,2));

cnzt=cnz;
[~,ddt1]=ddd1(z,r,phi,cnzt,cnr,cm);
cnzt=z+step*ddt1.cnz*beta(1,1);
[~,ddt2]=ddd1(z,r,phi,cnzt,cnr,cm);
cnzt=z+step*(ddt1.cnz*beta(1,2)+ddt2.cnz*beta(2,2));
[~,ddt3]=ddd1(z,r,phi,cnzt,cnr,cm);
cnzt=z+step*(ddt1.cnz*beta(1,3)+ddt2.cnz*beta(2,3)+ddt3.cnz*beta(3,3));
[~,ddt4]=ddd1(z,r,phi,cnzt,cnr,cm);
cnzt=z+step*(ddt1.cnz*beta(1,4)+ddt2.cnz*beta(2,4)+ddt3.cnz*beta(3,4)+ddt4.cnz*beta(4,4));
[~,ddt5]=ddd1(z,r,phi,cnzt,cnr,cm);
cnzt=z+step*(ddt1.cnz*beta(1,5)+ddt2.cnz*beta(2,5)+ddt3.cnz*beta(3,5)+ddt4.cnz*beta(4,5)+ddt5.cnz*beta(5,5));
[~,ddt6]=ddd1(z,r,phi,cnzt,cnr,cm);
ydelcnz=step*(ddt1.cnz*gama(1,2)+ddt2.cnz*gama(2,2)+ddt3.cnz*gama(3,2)+ddt4.cnz*gama(4,2)+ddt5.cnz*gama(5,2)+ddt6.cnz*gama(6,2));

cnrt=cnr;
[~,ddt1]=ddd1(z,r,phi,cnz,cnrt,cm);
cnrt=z+step*ddt1.cnr*beta(1,1);
[~,ddt2]=ddd1(z,r,phi,cnz,cnrt,cm);
cnrt=z+step*(ddt1.cnr*beta(1,2)+ddt2.cnr*beta(2,2));
[~,ddt3]=ddd1(z,r,phi,cnz,cnrt,cm);
cnrt=z+step*(ddt1.cnr*beta(1,3)+ddt2.cnr*beta(2,3)+ddt3.cnr*beta(3,3));
[~,ddt4]=ddd1(z,r,phi,cnz,cnrt,cm);
cnrt=z+step*(ddt1.cnr*beta(1,4)+ddt2.cnr*beta(2,4)+ddt3.cnr*beta(3,4)+ddt4.cnr*beta(4,4));
[~,ddt5]=ddd1(z,r,phi,cnz,cnrt,cm);
cnrt=z+step*(ddt1.cnr*beta(1,5)+ddt2.cnr*beta(2,5)+ddt3.cnr*beta(3,5)+ddt4.cnr*beta(4,5)+ddt5.cnr*beta(5,5));
[~,ddt6]=ddd1(z,r,phi,cnz,cnrt,cm);
ydelcnr=step*(ddt1.cnr*gama(1,2)+ddt2.cnr*gama(2,2)+ddt3.cnr*gama(3,2)+ddt4.cnr*gama(4,2)+ddt5.cnr*gama(5,2)+ddt6.cnr*gama(6,2));

cmt=cm;
[~,ddt1]=ddd1(z,r,phi,cnz,cnr,cmt);
cmt=z+step*ddt1.cm*beta(1,1);
[~,ddt2]=ddd1(z,r,phi,cnz,cnr,cmt);
cmt=z+step*(ddt1.cm*beta(1,2)+ddt2.cm*beta(2,2));
[~,ddt3]=ddd1(z,r,phi,cnz,cnr,cmt);
cmt=z+step*(ddt1.cm*beta(1,3)+ddt2.cm*beta(2,3)+ddt3.cm*beta(3,3));
[~,ddt4]=ddd1(z,r,phi,cnz,cnr,cmt);
cmt=z+step*(ddt1.cm*beta(1,4)+ddt2.cm*beta(2,4)+ddt3.cm*beta(3,4)+ddt4.cm*beta(4,4));
[~,ddt5]=ddd1(z,r,phi,cnz,cnr,cmt);
cmt=z+step*(ddt1.cm*beta(1,5)+ddt2.cm*beta(2,5)+ddt3.cm*beta(3,5)+ddt4.cm*beta(4,5)+ddt5.cm*beta(5,5));
[~,ddt6]=ddd1(z,r,phi,cnz,cnr,cmt);
ydelcm=step*(ddt1.cm*gama(1,2)+ddt2.cm*gama(2,2)+ddt3.cm*gama(3,2)+ddt4.cm*gama(4,2)+ddt5.cm*gama(5,2)+ddt6.cm*gama(6,2));

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
        
  