function [new_value]=calc_interp2(old_value)
intx1=1:129;
intx2=1:0.1:129;
inty1=1:129;
inty2=1:0.1:129;
[intx1,inty1]=meshgrid(intx1,inty1);
[intx2,inty2]=meshgrid(intx2,inty2);
new_value=interp2(intx1,inty1,old_value,intx2,inty2);