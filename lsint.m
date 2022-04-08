%ls=3/16*pi^0.5*df*(x0./k0)^2*(wpi^2./(w0*wce))^2*x0;
function [lsreci]=lsint(deltan,k0,wpe,w,wce,ksi0)
% syms x;
lsreci=3/16*pi^0.5*deltan^2*(ksi0/k0)^2*(wpe^2./(w *wce))^2*ksi0;   
% tao=int(3/16*pi^0.5*deltan^2*(x./k0)^2*(wpe^2./(w *wce))^2*x,x,0,x2); 