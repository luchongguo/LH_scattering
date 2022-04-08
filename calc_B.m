function [B]=calc_B(equil)
ii=120;
jj=30;
B=[equil.ptBx(ii,jj),equil.ptBy(ii,jj),equil.ptBPHI(ii,jj)];