function [z,r,phi,cnz,cnr,cm,us]=rk11(z,r,phi,cnz,cnr,cm,us)
[~,ddt,ii,jj]=ddd1(z,r,phi,cnz,cnr,cm);
step=5e-8;
z=z+ddt.z*step;
r=r+ddt.r*step;
phi=phi+ddt.phi*step;
cnz=ddt.cnz*step;
cnr=ddt.cnr*step;
cm=ddt.cm*step;
us=us+sqrt((ddt.z*step)^2+(ddt.r*step)^2+(ddt.phi*step)^2);
