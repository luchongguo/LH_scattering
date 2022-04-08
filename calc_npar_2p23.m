npar=-2.84;
ktheta11_1=[0.5 2 1]%input('Please input ktheta:');
%npar=input('Please input Npar:');
for i_kth=1:length(ktheta11_1)
    ktheta11=ktheta11_1(i_kth)

        calc_loop_old;
    save tmp.mat  i_kth ktheta11_1 npar;
    clear all;
    load tmp.mat
end