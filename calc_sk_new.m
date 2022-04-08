str = '/gpfs/scratch/chbwu/density_fluction_to_shenma/origindata/90331/simulate_new_20200409/';
files = dir(strcat(str,'*.mat'));
npar_1=[-2.23]
ktheta11=1;%input('Plse input ktheta:');
%npar=input('Please iut Npar:');

for i_npar=1:length(npar_1)
    

    
    for i_num=1:length(files)
        npar=npar_1(i_npar)
        load([str,files(i_num).name])
        calc_loop_old_3;
        save tmp.mat npar_1 i_npar ktheta11 str files i_num;
        clear all;
        load tmp.mat;

    end

end