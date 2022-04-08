function [new_value]=calc_interp(old_value)
int1=1:129;
int2=1:0.1:129;
new_value=interp1(int1,old_value,int2);