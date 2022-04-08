function [ii,jj]=findposition(z,r)
global eqdsk_r1 eqdsk_z1 
[~,ii]=min(abs(z-eqdsk_z1));
[~,jj]=min(abs(r-eqdsk_r1));
