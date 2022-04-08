for j=1:nray
    npar=-2.04;
    a=0.45;
    ksi0roui=0.1;
    cirho=1;

    z=0;
    r=2.3;
    phi=0;
    [imax,imin,jmax,jmin]=calc_position(z,r,eqdsk_z,eqdsk_r);
    calc_position_new;
%     [~,ii]=min(abs(z-eqdsk_z1));
%     [~,jj]=min(abs(r-eqdsk_r1));


%%
[xar1,xar2,yar1,yar2,npar_acc]=calc_xyar(Ne1(ii,jj),Ni1(ii,jj),Btot1(ii,jj));
[eps]=calc_eps1(xar1,xar2,yar1,yar2);
[nperps,nperpf,D]=calc_nperpsf(eps,npar);
%[eps,nperps(1),nperpf(1),D]=calc_eps(xar1(1),xar2(1),yar1(1),yar2(1),npar);
kpar=w*npar./c;

if imag(nperps(1))==0
    nperp=nperps;
end
 
% B_theta=(Bz.*dpsidr-Br.*dpsidz)./gradpsi;
cnteta1(1,j)=0;
cnphi=npar*Btot1(ii,jj)/Bt1(ii,jj);
cn2=npar^2+nperp^2;
cnrho2=cn2-cnteta1(1,j)^2-cnphi^2;
cnrho=sqrt(cnrho2);
[dpsidr1]=calc_interp22(imin,imax,jmin,jmax,dpsidr);
[dpsidz1]=calc_interp22(imin,imax,jmin,jmax,dpsidz);
gradpsi(ii,jj)=sqrt(dpsidr1(ii,jj).^2+dpsidz1(ii,jj).^2);
cnr=(cirho*cnrho*dpsidr1(ii,jj))/gradpsi(ii,jj);
cnz=(cirho*cnrho*dpsidz1(ii,jj))/gradpsi(ii,jj);
cm=cnphi*eqdskr(jj);

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
for i=1:lengthray
    complex(j,i)
   
%     [ddd,ddt,ii,jj]=ddd1(z,r,phi,cnz,cnr,cm);
%     dddd(i,j)=ddd;
%     ddt1(i,j)=ddt;
%     ii1(i,j)=ii;
%     jj1(i,j)=jj;
%     sne1(i,j)=Ne1(ii,jj);
%     Btot(i,j)=sqrt(Br1(ii,jj)^2+Bz1(ii,jj)^2+Bt1(ii,jj)^2);
    
%%
%if (Ne1(ii,jj)==0 || r<eqdsk_r(1) || r>eqdsk_r(end) || ddt.z>c || ddt.r>c) %ddt>c 属于什么问题？
%if inpolygon(r,z,gvar.rlim,gvar.zlim) && inpolygon(r,z,gvar.rbbbs,gvar.zbbbs)==0
    if inpolygon(r,z,gvar.rlim,gvar.zlim)
    [ddd,ddt,ii,jj]=calc_ddd1(z,r,phi,cnz,cnr,cm);
    if abs(ddt.z)<c & abs(ddt.r)<c
    dddd1(i,j)=ddd;
    ddt1(i,j)=ddt;
    ii1(i,j)=ii;
    jj1(i,j)=jj;
    sne1(i,j)=Ne1(ii,jj);
    Btot(i,j)=sqrt(Br1(ii,jj)^2+Bz1(ii,jj)^2+Bt1(ii,jj)^2);
    [xar1,xar2,yar1,yar2,npar_acc]=calc_xyar(sne1(i,j),sne1(i,j),Btot(i,j));
    [eps]=calc_eps1(xar1,xar2,yar1,yar2);
    [nperps,nperpf,D]=calc_nperpsf(eps,Npar1(i,j));
    npar_acc1(i,j)=npar_acc;
    sbr(i,j)=Br1(ii,jj);
    sbz(i,j)=Bz1(ii,jj);
    sbt(i,j)=Bt1(ii,jj);
    nperps1(i,j)=imag(nperps);
    nperps11(i,j)=nperps;
    nperpf1(i,j)=imag(nperpf);
    ste(i,j)=Te1(ii,jj);
    sti(i,j)=Te1(ii,jj);
%     [cnprim]=absorb(nperps1(i,j),Npar1(i,j),ste(i,j),sti(i,j),sne1(i,j),sne1(i,j),Btot(i,j),Zeff,w);
%     cnprim1.e(i,j)=cnprim.e; 
%     cnprim1.i(i,j)=cnprim.i;
%     cnprim1.cl(i,j)=cnprim.cl;
%     deltaP.e(i,j)=exp(-2*cnprim1.e(i,j)*tstep);
%     deltaP.i(i,j)=exp(-2*cnprim1.i(i,j)*tstep);
%     deltaP.cl(i,j)=exp(-2*cnprim1.cl(i,j)*tstep);

%%
[z,r,phi,cnz,cnr,cm,us]=rk2(z,r,phi,cnz,cnr,cm,us);
%[ii,jj]=findposition(z,r);
% [~,ii]=min(abs(z-eqdsk_z1));
% [~,jj]=min(abs(r-eqdsk_r1));
[imax,imin,jmax,jmin]=calc_position(z,r,eqdsk_z,eqdsk_r);
calc_position_new;
z1(i+1,j)=z;
r1(i+1,j)=r;
phi1(i+1,j)=phi;
cnz1(i+1,j)=cnz;
cnr1(i+1,j)=cnr;
cm1(i+1,j)=cm;
us1(i+1,j)=us;
[teta]=calc_teta11(z1(i+1,j),r1(i+1,j),z1(i,j),r1(i,j),gvar);
teta1(i+1,j)=teta1(i,j)+teta;
% if teta1(i,j)>teta
%     teta1(i+1,j)=teta+2*pi;
% else
%     teta1(i+1,j)=teta;
% end

Npar1(i+1,j)=(Br1(ii,jj)*cnr+Bz1(ii,jj)*cnz+Bt1(ii,jj)*cm/eqdskr(jj))/Btot(i,j);
%


%%  区分刮削层和内部涨落水平
if  j~=1
    Deltan11(i,j)=Deltan1(ii,jj);
    if  Deltan1(ii,jj)~=0  &&  (imag(nperps)~=0 && imag(nperpf)~=0)~=1
        [snperm,ds,cnz,cnr,cm,xrandom]=densfluc(z,r,phi,cnz,cnr,cm,tstep);       %在这里并没有给出新的（cnz,cnr,cm）,已修改
        xrandom1(i,j)=xrandom;
        snperm1(i,j)=snperm;
       
    end
% else
%     ds=0;
end
%  us=us+ds;      %意义在哪？如果(z,r)不变的话，那么ds的变化是怎么来的？

%% 
%[z,r,phi,cnz,cnr,cm,us]=rk11(z,r,phi,cnz,cnr,cm,us);
%[z,r,phi,cnz,cnr,cm,us]=rk2(z,r,phi,cnz,cnr,cm,us);



%     if ((abs(Npar1(i,j))<npar_acc1(i,j))||(imag(nperps)>0&&imag(nperpf)>0))     %截止条件
% %         Npar1=Npar(i,j);
% %         z1(i,j)=z;
% %         r1(i,j)=r;
% %         phi1(i,j)=phi;
% %         us1(i,j)=us;
% %         cnz1(i,j)=cnz;
% %         cnr1(i,j)=cnr;
% %         cm1(i,j)=cm;
% %         rhopsi1(i,j)=rhopsi11(ii,jj);
% %         cnteta(i,j)=(cnr*Br1(ii,jj)+cnz*Bz1(ii,jj))/sqrt(Br1(ii,jj)^2+Bz1(ii,jj)^2);
%         rhopsi1(i+1,j)=rhopsi11(ii,jj);
%         cnteta1(i+1,j)=(cnr*Br1(ii,jj)+cnz*Bz1(ii,jj))/sqrt(Br1(ii,jj)^2+Bz1(ii,jj)^2);
%         break;
%     else
%         if rhopsi11(ii,jj)>1.2      %在rho为1.2处发生反射
%             cnteta=(cnr*Br1(ii,jj)+cnz*Bz1(ii,jj))/sqrt(Br1(ii,jj)^2+Bz1(ii,jj)^2);
%             cnpsi=(cnz*Br1(ii,jj)-cnr*Bz1(ii,jj))/sqrt(Br1(ii,jj)^2+Bz1(ii,jj)^2);
%             cnzref=(cnteta*Bz1(ii,jj)-cnpsi*Br1(ii,jj))/sqrt(Br1(ii,jj)^2+Bz1(ii,jj)^2);
%             cnrref=(cnteta*Br1(ii,jj)+cnpsi*Bz1(ii,jj))/sqrt(Br1(ii,jj)^2+Bz1(ii,jj)^2);
%             cnz=cnzref;
%             cnr=cnrref;
%         end
%         rhopsi1(i+1,j)=rhopsi11(ii,jj);
%         cnteta1(i+1,j)=(cnr*Br1(ii,jj)+cnz*Bz1(ii,jj))/sqrt(Br1(ii,jj)^2+Bz1(ii,jj)^2);
%     end
    end
end
%%
end
end