%testDir='F:\cbwu\锟届公\density_fluction_0408\origindata\';
shotno=90331
testDir=['E:\OneDrive - mail.ustc.edu.cn\个人信息\学籍相关\个人论文\密度涨落\LH_scattering_code\origindata\',num2str(shotno),'\'];

%gfile = fullfile(testDir,'g080307.040000_old');
%gfile = fullfile(testDir,'g071320.02800');

lam_ne=0.1;
lam_te=0.05;
lam_delta=0.2;

gfile = fullfile(testDir,'g090331.006110');
load([testDir,'profiles.mat']);
load([testDir,'nedr_106818.mat']);
load([testDir,'Bdata.mat']);
read_gfile;
[eqdsk_r1]=calc_interp(eqdsk_r);
[eqdsk_z1]=calc_interp(eqdsk_z);
[rhopsi11]=calc_interp2(rhopsi);

[ne1,te1,delta1,rho1,ktheta1]=calc_edge_20211230(nedr,eqdsk_r1,rhopsi11,ne,te,lam_ne,lam_te,lam_delta);
[Ne]=calc_profiles(ne,ne1,rho,rhopsi11,rho1);
[Te]=calc_profiles(te,te1,rho,rhopsi11,rho1);
[Ti]=calc_profiles(te/TeTi,te1/TeTi,rho,rhopsi11,rho1);
Delta_11=[zeros(1,41) 5.575*(rho(42:51)-0.8).^2];
[Deltan]=calc_profiles(Delta_11,delta1,rho,rhopsi11,rho1);
[Ktheta]=calc_profiles(zeros(1,51),ktheta1,rho,rhopsi11,rho1);

figure
plot(rhopsi11(641,641:end),Ne(641,641:end))
figure
plot(rhopsi11(641,641:end),Deltan(641,641:end))
figure
[rgrid1]=calc_interp(rgrid);
plot(rgrid1(1000:end)*1000,Ne(641,1000:end)./10e17)

