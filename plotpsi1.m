function []=plotpsi1(gvar,prof)


%%
X = linspace(0, gvar.rdim, gvar.nw);
Z = linspace(0, gvar.zdim, gvar.nw);
rgrid = X+gvar.rleft;
zgrid = Z-(gvar.zmid+gvar.zdim/2);
psirz = gvar.psirz;
ssimag = gvar.simag;    %磁轴处psi值
ssibry = gvar.sibry;    %    最外闭合磁面处psi值
pres = gvar.pres;
%%
a=(max(gvar.rbbbs)-min(gvar.rbbbs))/2;
b=(max(gvar.zbbbs)-min(gvar.zbbbs))/2;
R=(max(gvar.rbbbs)+min(gvar.rbbbs))/2;
k=b/a;
dtriu=(R-gvar.rbbbs(find(gvar.zbbbs==max(gvar.zbbbs))))/a;
dtrid=(R-gvar.rbbbs(find(gvar.zbbbs==min(gvar.zbbbs))))/a;
rhopsi=(psirz'-ssimag)/(ssibry-ssimag);

% rhopsip=interp2(rgrid,zgrid,rhopsi,rgrid,zgrid,'cubic');


%%
figure(1);
contourf(rgrid,zgrid,prof,50);
% %  figure();
% %  [c1,h1]=contourf(rgrid,zgrid,rhopsi,50);
% [r,z]=meshgrid(rgrid,zgrid);
%[c,h]=contour(X,Y,Z,'LevelList',[5.0]);
%contourf(rgrid,zgrid,Bp,50);
hold on;
plot(gvar.rbbbs,gvar.zbbbs,'r','LineWidth',2);   %最外闭合磁面
plot(gvar.rlim,gvar.zlim,'b','LineWidth',2);    %限制器
colorbar();
axis equal;
xlabel('R(m)');
ylabel('Z(m)');
hold on;