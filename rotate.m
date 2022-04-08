function [cnzm,cnrm,cmm]=rotate(cnz,cnr,cm,r,br,bz,bphi,cnpar,skperp,skperm,w,beta1)
%(cnz,cnr,cm)决定了cnpar，进而决定了snperp，snperm，故会存在一些偏差，真正要做到beta=0时，skz=skzm，要满足的条件比较多

c=299792458;
btot=sqrt(br^2+bz^2+bphi^2);
%skpr=(cnz*bz+cnr*br+cm/r*bphi)/btot;
skpr=cnpar*w/c;
skz=cnz*w/c;
skr=cnr*w/c;
skph=cm/r*w/c;

pt1=skph-skpr*bphi/btot;
pt2=(skz*br-skr*bz)/btot;
pt3=(skph*br-skr*bphi)/btot;
pt4=skperp*skperm*cos(beta1)+skpr^2;
skphm=skpr*bphi/btot+pt1*cos(beta1)*skperm/skperp+pt2*sin(beta1)*skperm/skperp;
skzm=(br*pt4/btot-skr*skpr-skphm*pt3)/pt2;
skrm=btot*(skpr-skzm*bz/btot-skphm*bphi/btot)/br;
cnzm=c/w*skzm;
cnrm=c/w*skrm;
cmm=r*c/w*skphm;







