function [nperps,nperpf,D]=calc_nperpsf(eps,npar)



A=eps.perp;
B=(-eps.perp+npar^2)*(eps.perp+eps.par)+eps.xy^2;
C=eps.par*((eps.perp-npar^2)^2-eps.xy^2);

nperps2=(-B+(B^2-4*eps.perp*C)^0.5)/2/A;
nperpf2=(-B-(B^2-4*eps.perp*C)^0.5)/2/A;

nperps=sqrt(nperps2);
nperpf=sqrt(nperpf2);
D=A*nperps2^2+B*nperps2+C;

