function [depsdxy]=calc_depsdxy(rho,xar,yar)
% depsdxy.perpx=1/yar(:)^2;
% depsdxy.parx=-1;
% depsdxy.xyx=1/yar(:);
% depsdxy.perpy=-2*xar(:)/yar(:)^3;
% depsdxy.pary=0;
% depsdxy.xyy
step=1e-5;

[epsp]=calc_eps(xar,yar);
rhom=rho+step;
[epsm]=calc_eps(xar,yar);

