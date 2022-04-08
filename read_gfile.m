gfileVar = readGfile(gfile); 
gvar = gfileVar;
X = linspace(0, gvar.rdim, gvar.nw);
Z = linspace(0, gvar.zdim, gvar.nw);
rgrid = X+gvar.rleft;
zgrid = Z-(gvar.zmid+gvar.zdim/2);
% figure(1);
% contourf(rgrid,zgrid,-gvar.psirz',50);
% hold on;
% plot(gvar.rbbbs,gvar.zbbbs,'r','LineWidth',2);   %最外闭合磁面
% plot(gvar.rlim,gvar.zlim,'b','LineWidth',2);    %限制器
% colorbar();
% axis equal;
% xlabel('R(m)');
% ylabel('Z(m)');
[dpsidr,dpsidz]=gradient(gvar.psirz');
X = linspace(0, gvar.rdim, gvar.nw);
Z = linspace(0, gvar.zdim, gvar.nw);
eqdsk_r = X+gvar.rleft;
eqdsk_z = Z-(gvar.zmid+gvar.zdim/2);
% hold on;
 rhopsi=(gvar.psirz'-gvar.simag)/(gvar.sibry-gvar.simag);
% 
% contourf(rgrid,zgrid,rhopsi,50);