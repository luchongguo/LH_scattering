clear;clc;
format long
data_basic;
%shotno=input('Please input the number of ray:')
shotno=90331;
w11=input('Please input wave frequency:');
w1=w11*2*pi*10^9; %���Ӳ���ԴƵ��
data_prepare;
data_global_old;
tstep=1e-11;


%tic
%%
% data_input;
% calc_loop_old;
% data_global; 
% calc_loop;
%%
