[eqdskz]=calc_interp1(imax,imin,eqdsk_z);
[eqdskr]=calc_interp1(jmax,jmin,eqdsk_r);
[Bt1]=calc_interp22(imin,imax,jmin,jmax,Bt);
[Bp1]=calc_interp22(imin,imax,jmin,jmax,Bp);
[Br1]=calc_interp22(imin,imax,jmin,jmax,Br);
[Bz1]=calc_interp22(imin,imax,jmin,jmax,Bz);
[Ne1]=calc_interp22(imin,imax,jmin,jmax,Ne);
[Ne1]=calc_interp22(imin,imax,jmin,jmax,Ne);
[Te1]=calc_interp22(imin,imax,jmin,jmax,Te);
Ni1=Ne1/1.25;
[~,ii]=min(abs(r-eqdskr));
[~,jj]=min(abs(z-eqdskz));
Btot1(ii,jj)=sqrt(Br1(ii,jj)^2+Bz1(ii,jj)^2+Bt1(ii,jj)^2);

[Deltan1]=calc_interp22(imin,imax,jmin,jmax,Deltan);
[rhopsi11]=calc_interp22(imin,imax,jmin,jmax,rhopsi);