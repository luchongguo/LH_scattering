function [rho,teta]=calc_teta(z,r,gvar)
zpos=z-gvar.zmaxis;
rpos=r-gvar.rmaxis;
pos=sqrt(zpos.^2+rpos.^2);
rho=pos./sqrt(zpos(end)^2+rpos(end)^2);      %²»Ì«ÕýÈ·
%rrho=rpos./rpos(end);
% if zpos>0
%     teta=acos(rpos./pos);
% else
%     teta=-acos(rpos./pos);
% end
if zpos>0 && rpos>0
    teta=acos(rpos./pos);
elseif zpos>0 &&rpos<0
    teta=acos(rpos./pos);
elseif zpos<0 &&rpos>0
    teta=2*pi-acos(rpos./pos);
else
    teta=2*pi-acos(rpos./pos);
end