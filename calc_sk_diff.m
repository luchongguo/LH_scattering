str = '/gpfs/scratch/chbwu/density_fluction_to_shenma/origindata/90331/simulate_new_20200409/';
files = dir(strcat(str,'*.mat'));
npar_1=[-2.23]
ii_ktheta11=[0.5 2];%input('Plse input ktheta:');
%npar=input('Please iut Npar:');
for i_npar=1:length(npar_1)
    

    
    for i_num=1:length(ii_ktheta11)
        npar=npar_1(i_npar)
        ktheta11=ii_ktheta11(i_num)
        load([str,files(4).name])
        calc_loop_old_3;
        save tmp.mat npar_1 i_npar ktheta11 str files i_num ii_ktheta11;
        clear all;
        load tmp.mat;

    end

end