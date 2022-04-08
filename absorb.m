function [cnprim]=absorb(cnper,cnpar,Te,Ti,ne,ni,Btot,Zeff,w,ddd)
%ve 要注意单位

data_basic;  
Ti=Ti*1000*1.6e-19;
Te=Te*1000*1.6e-19;
[xar1,xar2,yar1,yar2,~]=calc_xyar(ne,ni,Btot);
wpe=(ne*e*e/ee/m).^0.5; 
ve=(sqrt(Te/m));
vi=(sqrt(Ti/M));
xoe=c/(sqrt(2)*cnpar*ve);   
xoi=c/(sqrt(2)*cnpar*vi);
di.e=2*pi^0.5*xar1*cnpar^2*cnper^2*xoe^3*exp(-xoe^2);
di.i=2*pi^0.5*xar2*w*cnper^4*xoi^3*exp(-xoi^2);
nu0=wpe^4*16/2/pi/ne/ve^3;
nuei=2/3*sqrt(pi)*nu0*Zeff;
di.cl=nuei/w*(xar1/yar1^2*cnper^2+xar1*cnpar^2)*cnper^2;


cnpars=cnpar^2;
cnper2=cnper^2;
cnper4=cnper^4;

[eps]=calc_eps1(xar1,xar2,yar1,yar2);


p0=eps.par*((cnpars-eps.perp)^2-eps.xy^2);
p2=(eps.perp+eps.par)*(cnpars-eps.perp)+eps.xy^2;
p4=eps.perp;
p6=-3/2*xar2*(vi/c)^2-3/8*xar1/yar1^4*(ve/c);

% dddcnper=(6.d0*p6*cnper4+4*p4*cnper2+2*p2)*cnper;
% cnprim.e=-(di.e/dddcnper);
% cnprim.i=-(di.i/dddcnper);
% cnprim.cl=-(di.cl/dddcnper);
%dddcnper=(6.d0*p6*cnper4+4*p4*cnper2+2*p2)*cnper;
cnprim.e=-(di.e/ddd.w);
cnprim.i=-(di.i/ddd.w);
cnprim.cl=-(di.cl/ddd.w);
cnprim.e=abs(cnprim.e);
cnprim.i=abs(cnprim.i);
cnprim.cl=abs(cnprim.cl);
