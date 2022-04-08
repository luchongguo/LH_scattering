function [D0,P0246]=calc_P0246(eps,npar,nperp)

%%
npars=npar*npar;
nperp2=nperp*nperp;
nperp4=nperp2*nperp2;
nperp6=nperp2*nperp4;

P0246.p0=eps.par*((npars-eps.perp)^2-eps.xy^2);
P0246.p2=(eps.perp+eps.par)*(npars-eps.perp)+eps.xy^2;
P0246.p4=eps.perp;
%P0246.p6=-3/2*xar2*(vi/c)^2-3/8*xar1/yar1^4*(ve/c);
P0246.p6=0;
D0=P0246.p6*nperp6+P0246.p4*nperp4+P0246.p2*nperp2+P0246.p0;
