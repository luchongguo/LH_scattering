function [new_value]=calc_interp1(imax,imin,old_value)
int1=imin:imax;
int2=imin:0.01:imax;
if imin~=imax
    new_value=interp1(int1,old_value(imin:imax),int2);
else
    new_value=old_value(imax);
end