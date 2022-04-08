function [new_value]=calc_interp22(imin,imax,jmin,jmax,old_value)

if imin~=imax & jmin~=jmax
    intx1=imin:imax;
    intx2=imin:0.01:imax;
    inty1=jmin:jmax;
    inty2=jmin:0.01:jmax;
    [intx1,inty1]=meshgrid(intx1,inty1);
    [intx2,inty2]=meshgrid(intx2,inty2);
    new_value=interp2(intx1,inty1,old_value(imin:imax,jmin:jmax),intx2,inty2);

elseif imin==imax & jmin~=jmax
    intj1=jmin:jmax;
    intj2=jmin:0.01:jmax;
    new_value=interp1(intj1,old_value(imin:imax,jmin:jmax),intj2);
    new_value=repmat(new_value,length(intj2),1);

elseif jmin==jmin & imin~=imax
    inti1=imin:imax;  
    inti2=imin:0.01:imax;
    new_value=interp1(inti1,old_value(imin:imax,jmin:jmax),inti2);
    new_value=repmat(new_value,1,length(inti2));

elseif imin==imax & jmin==jmax
    intx2=imin:0.01:imax;
    inty2=jmin:0.01:jmax;
    new_value=repmat(old_value(imin:imax,jmin:jmax),length(intx2),length(inty2));
 end

    
    
    