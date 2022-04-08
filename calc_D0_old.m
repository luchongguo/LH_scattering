function [D0,cnperp]=calc_D0_old(xar1,xar2,yar1,yar2,gamma1,cn,ii,jj)

global Te1 Ne1 Btot1 Ti1
M=2*1.67*10^-27;
m=9.1*10^-31;
c=299792458;
%%
Ti=Ti1(ii,jj)*1000*1.6e-19;
Te=Te1(ii,jj)*1000*1.6e-19;
[eps]=calc_eps1(xar1,xar2,yar1,yar2);

dc=cos(gamma1);
ds=sin(gamma1);
cnpar=cn*dc;
cnperp=cn*ds;

vi=sqrt(2*Ti/M);
ve=sqrt(2*Te/m);

%%
cnpars=cnpar*cnpar;
cnperp2=cnperp*cnperp;
cnperp4=cnperp2*cnperp2;
cnperp6=cnperp2*cnperp4;

P0246.p0=eps.par*((cnpars-eps.perp)^2-eps.xy^2);
P0246.p2=(eps.perp+eps.par)*(cnpars-eps.perp)+eps.xy^2;
P0246.p4=eps.perp;
P0246.p6=-3/2*xar2*(vi/c)^2-3/8*xar1/yar1^4*(ve/c);
%P0246.p6=0;
D0=P0246.p6*cnperp6+P0246.p4*cnperp4+P0246.p2*cnperp2+P0246.p0;
%%






