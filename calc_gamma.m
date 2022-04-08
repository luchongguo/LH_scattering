function [gamm,cn]=calc_gamma(bz,br,bphi,cnz,cnr,cm,r)
cn=sqrt(cnz^2+cnr^2+cm^2/r^2);
gg=cnz*bz+cnr*br+bphi*cm/r;
btot=sqrt(bphi^2+bz^2+br^2);
arg=gg/(btot*cn);
if arg>1|arg==1
    arg=1
    gamm=0
elseif arg<-1 |arg==-1
    arg=-1;
    gamm=pi;
else
    gamm=acos(arg);
end