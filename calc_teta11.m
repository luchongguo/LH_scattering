function [teta]=calc_teta11(zstart,rstart,zend,rend,gvar)
zzstart=zstart-gvar.zmaxis;
rrstart=rstart-gvar.rmaxis;
zzend=zend-gvar.zmaxis;
rrend=rend-gvar.rmaxis;
zz=zend-zstart;
rr=rend-rstart;
l1=sqrt(zzstart^2+rrstart^2);
l2=sqrt(zzend^2+rrend^2);
l3=sqrt(zz^2-rr^2);
costeta=(l1^2+l2^2-l3^2)/(2*l1*l2);
teta=acos(costeta);
teta=abs(teta);


%rrho=rpos./rpos(end);
% if zpos>0
%     teta=acos(rpos./pos);
% else
%     teta=-acos(rpos./pos);
% end