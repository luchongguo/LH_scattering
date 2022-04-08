npar=[-2.84];
ktheta11=2;%input('Please input ktheta:');
%npar=input('Please input Npar:');
Kmean=[-0.5 0.5]
for i_npar=1:length(Kmean)
    kmean=Kmean(i_npar)

        calc_loop_old_3;
    save tmp.mat Kmean i_npar ktheta11 npar;
    clear all;
    load tmp.mat
end