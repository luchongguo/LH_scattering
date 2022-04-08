% global Bt1 Bp1 Br1 Bz1 eqdsk_r1 eqdsk_z1  Ne1 Ni1 Te1 w Btot1 Deltan1
% clear global Bt1;
% clear global Bp1;
% clear global Br1;
% clear global Bz1;
% clear global eqdsk_r1;
% clear global eqdsk_z1;
% clear global Ne1;
% clear global Ni1;
% clear global Te1;
% clear global w;
% clear global Btot1;
% clear global Deltan1;
clear global;

%%
%global Bt Bp Br Bz eqdsk_r eqdsk_z  Ne Ni Te w  Deltan rhopsi 
global Bt1 Bp1 Br1 Bz1 eqdsk_r1 eqdsk_z1  Ne1 Ni1 Te1 Ti1 w Btot1 Deltan1 Ktheta1
[Bt1]=calc_interp2(Bt);
[Bp1]=calc_interp2(Bp);
[Br1]=calc_interp2(Br);
[Bz1]=calc_interp2(Bz);
[eqdsk_r1]=calc_interp(eqdsk_r);
[eqdsk_z1]=calc_interp(eqdsk_z);
[Ne1]=Ne;
[Te1]=Te;
[Ti1]=Ti;
[Deltan1]=Deltan;
[Ktheta1]=Ktheta;
Btot1=(Br1.^2+Bz1.^2+Bt1.^2).^0.5;
Ni1=Ne;
%[rhopsi11]=rhopsi;
w=w1;
[dpsidr1]=calc_interp2(dpsidr);
[dpsidz1]=calc_interp2(dpsidz); 

