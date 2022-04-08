%function [rho,teta]=calc_teta1(eqdsk_r,eqdsk_z,gvar,a,k)
[r,z]=meshgrat(eqdsk_r,eqdsk_z);
zpos=z-gvar.zmaxis;

rpos=r-gvar.rmaxis;
pos=sqrt(zpos.^2/(a*k)^2+rpos.^2/a^2);
%rho=pos/sqrt(zpos(end)^2+rpos(end)^2);      %不太正确
%rho=pos/sqrt(zpos(end)^2+rpos(end)^2);      %不太正确
%rrho=rpos./rpos(end);
% if zpos>0
%     teta=acos(rpos./pos);
% else
%     teta=-acos(rpos./pos);
% end
for i=1:129
    for j=1:129
        if zpos(i,j)>0 && rpos(i,j)>0
            teta(i,j)=acos(rpos(i,j)/pos(i,j));
        elseif zpos(i,j)>0 &&rpos(i,j)<0
            teta(i,j)=acos(rpos(i,j)/pos(i,j));
        elseif zpos(i,j)<0 &&rpos(i,j)>0
            teta(i,j)=2*pi-acos(rpos(i,j)/pos(i,j));
        else
            teta(i,j)=2*pi-acos(rpos(i,j)/pos(i,j));
        end
    end
end
z1=pos.*sin(teta);
r1=pos.*cos(teta);
z1=z1+gvar.zmaxis;
r1=r1+gvar.rmaxis;