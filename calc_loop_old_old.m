data_zero;
for j=1:nray
   % npar=-2.04;
    a=0.45;
    ksi0roui=0.1;
    cirho=1;

    z=0;
    r=2.321;
    phi=0;
%     [imax,imin,jmax,jmin]=calc_position(z,r,eqdsk_z,eqdsk_r);
%     calc_position_new;
    [~,ii]=min(abs(z-eqdsk_z1));
    [~,jj]=min(abs(r-eqdsk_r1));


%%
[xar1,xar2,yar1,yar2,npar_acc]=calc_xyar(Ne1(ii,jj),Ni1(ii,jj),Btot1(ii,jj));
[eps]=calc_eps1(xar1,xar2,yar1,yar2);
[nperps,nperpf,D]=calc_nperpsf(eps,npar);
kpar=w*npar./c;

if imag(nperps(1))==0
    nperp=nperps;
end
 
cnteta1(1,j)=0;
cnphi=npar*Btot1(ii,jj)/Bt1(ii,jj);
cn2=npar^2+nperp^2;
cnrho2=cn2-cnteta1(1,j)^2-cnphi^2;
cnrho=sqrt(cnrho2);
gradpsi(ii,jj)=sqrt(dpsidr1(ii,jj).^2+dpsidz1(ii,jj).^2);
cnr=(cirho*cnrho*dpsidr1(ii,jj))/gradpsi(ii,jj);
cnz=(cirho*cnrho*dpsidz1(ii,jj))/gradpsi(ii,jj);
cm=cnphi*eqdsk_r1(jj);

npar1=(cnr*Br1(ii,jj)+cnz*Bz1(ii,jj)+Bt1(ii,jj)*cnphi)/Btot1(ii,jj);
us=0;
z1(1,j)=z;
r1(1,j)=r;
phi1(1,j)=phi;
cnz1(1,j)=cnz;
cnr1(1,j)=cnr;
cm1(1,j)=cm;
us1(1,j)=us;
Npar1(1,j)=npar1;
teta1(1,j)=0;



%%
for i=2:lengthray
    complex(j,i)
    
%%
%if (Ne1(ii,jj)==0 || r<eqdsk_r(1) || r>eqdsk_r(end) || ddt.z>c || ddt.r>c) %ddt>c 属于什么问题？
if inpolygon(r,z,gvar.rlim,gvar.zlim) && inpolygon(r,z,gvar.rbbbs,gvar.zbbbs)==0
  %  if inpolygon(r,z,gvar.rlim,gvar.zlim)
    [ddd,ddt,ii,jj]=calc_ddd1_old(z,r,phi,cnz,cnr,cm);
     if abs(ddt.z)<c & abs(ddt.r)<c
%     dddd1(i,j)=ddd;
%     ddt1(i-1,j)=ddt;
    ii1(i-1,j)=ii;
    jj1(i-1,j)=jj;
    ssne=Ne1(ii,jj);
    sne1(i-1,j)=Ne1(ii,jj);
    Btot(i-1,j)=sqrt(Br1(ii,jj)^2+Bz1(ii,jj)^2+Bt1(ii,jj)^2);
    BBtot=sqrt(Br1(ii,jj)^2+Bz1(ii,jj)^2+Bt1(ii,jj)^2);
    [xar1,xar2,yar1,yar2,npar_acc]=calc_xyar(ssne,ssne,BBtot);
    [eps]=calc_eps1(xar1,xar2,yar1,yar2);
    Npartmp=Npar1(i-1,j);
    [nperps,nperpf,D]=calc_nperpsf(eps,Npartmp);
    npar_acc1(i-1,j)=npar_acc;
    sbr(i-1,j)=Br1(ii,jj);
    sbz(i-1,j)=Bz1(ii,jj);
    sbt(i-1,j)=Bt1(ii,jj);
    nperps1(i-1,j)=imag(nperps);
    nperps11(i-1,j)=nperps;
    nperpf1(i-1,j)=imag(nperpf);
    nperpf11(i-1,j)=nperpf;
    ste(i-1,j)=Te1(ii,jj);
    stetmp=ste(i-1,j);
    sti(i-1,j)=Ti1(ii,jj);
    stitmp=sti(i-1,j);
[cnprim]=absorb(nperps11(i-1,j),Npar1(i-1,j),ste(i-1,j),sti(i-1,j),sne1(i-1,j),sne1(i-1,j),Btot(i-1,j),Zeff,w,ddd);
    cnprim1.e(i-1,j)=cnprim.e; 
    cnprim1.i(i-1,j)=cnprim.i;
    cnprim1.cl(i-1,j)=cnprim.cl;
%%
[z,r,phi,cnz,cnr,cm,us]=rk2_old(z,r,phi,cnz,cnr,cm,us);
%[ii,jj]=findposition(z,r);
[~,ii]=min(abs(z-eqdsk_z1));
[~,jj]=min(abs(r-eqdsk_r1));
% [imax,imin,jmax,jmin]=calc_position(z,r,eqdsk_z,eqdsk_r);
% calc_position_new;
z1(i,j)=z;
r1(i,j)=r;
phi1(i,j)=phi;
cnz1(i,j)=cnz;
cnr1(i,j)=cnr;
cm1(i,j)=cm;
us1(i,j)=us;
[teta]=calc_teta11(z1(i,j),r1(i,j),z1(i-1,j),r1(i-1,j),gvar);
teta1(i,j)=teta1(i,j)+teta;
Npar1(i,j)=(Br1(ii,jj)*cnr+Bz1(ii,jj)*cnz+Bt1(ii,jj)*cm/eqdsk_r1(jj))/Btot1(ii,jj);
%


%%  区分刮削层和内部涨落水平
if  j~=1
    Deltan11(i,j)=Deltan1(ii,jj);
    if  Deltan1(ii,jj)~=0  &&  (imag(nperps)~=0 && imag(nperpf)~=0)~=1
        [snperm,ds,cnz,cnr,cm,lent,lenl,beta1,Vgperm,nu1]=densfluc_old(z,r,phi,cnz,cnr,cm,tstep,ktheta11);       %在这里并没有给出新的（cnz,cnr,cm）,已修改
%         xrandom1(i,j)=xrandom;
         lent1(i,j)=lent;
         lenl1(i,j)=lenl;
         beta11(i,j)=beta1;
         Vgperm1(i,j)=Vgperm;
         nu11(i,j)=nu1;
         snperm1(i,j)=snperm;
       
    end
% else
%     ds=0;
end
%  us=us+ds;      %意义在哪？如果(z,r)不变的话，那么ds的变化是怎么来的？

%% 
%[z,r,phi,cnz,cnr,cm,us]=rk11(z,r,phi,cnz,cnr,cm,us);
%[z,r,phi,cnz,cnr,cm,us]=rk2(z,r,phi,cnz,cnr,cm,us);



    if ((abs(Npar1(i,j))<npar_acc1(i,j))||(imag(nperps)>0&&imag(nperpf)>0))     %截止条件

        rhopsi1(i,j)=rhopsi11(ii,jj);
        cnteta1(i,j)=(cnr*Br1(ii,jj)+cnz*Bz1(ii,jj))/sqrt(Br1(ii,jj)^2+Bz1(ii,jj)^2);
%         break;
    else
        if rhopsi11(ii,jj)>1.2      %在rho为1.2处发生反射
            cnteta=(cnr*Br1(ii,jj)+cnz*Bz1(ii,jj))/sqrt(Br1(ii,jj)^2+Bz1(ii,jj)^2);
            cnpsi=(cnz*Br1(ii,jj)-cnr*Bz1(ii,jj))/sqrt(Br1(ii,jj)^2+Bz1(ii,jj)^2);
            cnzref=(cnteta*Bz1(ii,jj)-cnpsi*Br1(ii,jj))/sqrt(Br1(ii,jj)^2+Bz1(ii,jj)^2);
            cnrref=(cnteta*Br1(ii,jj)+cnpsi*Bz1(ii,jj))/sqrt(Br1(ii,jj)^2+Bz1(ii,jj)^2);
            cnz=cnzref;
            cnr=cnrref;
        end
        rhopsi1(i,j)=rhopsi11(ii,jj);
        cnteta1(i,j)=(cnr*Br1(ii,jj)+cnz*Bz1(ii,jj))/sqrt(Br1(ii,jj)^2+Bz1(ii,jj)^2);
    end
    end
end
%%
end
end