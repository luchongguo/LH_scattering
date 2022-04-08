function [eps]=calc_eps1(xar1,xar2,yar1,yar2)
%%
%¹Ì¶¨ÊýÖµ
data_basic;  
%%
eps.perp=[];
eps.par=[];
eps.xy=[];
% eps.perp=1+xar1/(yar1)^2-xar2;
% eps.par=1-xar1-xar2;
% eps.xy=xar1/yar1;
eps.perp=1-xar1/(1-(yar1)^2)-xar2/(1-(yar2)^2);
eps.par=1-xar1-xar2;
eps.xy=xar1*yar1/(1-(yar1)^2);
%reps=[eps.perp,1i*eps.xy,0;-1i*eps.xy,eps.perp,0;0,0,eps.par];
