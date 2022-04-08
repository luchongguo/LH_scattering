function [imax,imin,jmax,jmin]=calc_position(z,r,eqdsk_z,eqdsk_r)

    z1=eqdsk_z-z;
    r1=eqdsk_r-r;
    for i=1:128
        if z1(i)*z1(i+1)<0
            imax=i+1;
            imin=i;
        elseif z1(i)==0
            imin=i;
            imax=i;
        end
    end
    for j=1:128
        if r1(j)*r1(j+1)<0
            jmax=j+1;
            jmin=j;
        elseif r1(j)==0
            jmin=j;
            jmax=j;
        end
    end

