%npar_1=[-1.82 -1.93 -2.1 -2.23 -2.37 -2.6 -2.84 -1.82 -1.93 -2.04 -2.26 -2.48]
npar_1=[-2.26]
ktheta11=1;%input('Plse input ktheta:');
%npar=input('Please iut Npar:');
for i_npar=1:length(npar_1)
    npar=npar_1(i_npar)

        calc_loop_old_3;
    save tmp.mat npar_1 i_npar ktheta11;
    clear all;
    load tmp.mat
end