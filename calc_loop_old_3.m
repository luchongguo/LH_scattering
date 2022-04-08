

tic;
npar=2.36
nray=1;
lengthray=100;
kmean=0;
tstep=2e-11;

load('E:\OneDrive - mail.ustc.edu.cn\¸öÈËÐÅÏ¢\Ñ§¼®Ïà¹Ø\¸öÈËÂÛÎÄ\ÃÜ¶ÈÕÇÂä\LH_scattering_code\origindata\90331\simulate_new_20200823_rho_0p8/data_f_2p45_90331_8e18_delta_0p2_new_inside_LCFS_20211230.mat')



    
for j=1:nray
j
   
    a=0.45;
    ksi0roui=0.1;
    cirho=-1;

    zz=0;

    rr=2.38;
    z=zz;
    r=rr;

    phi=0;

    [~,ii]=min(abs(z-eqdsk_z1));
    [~,jj]=min(abs(r-eqdsk_r1));

    [xar1,xar2,yar1,yar2,npar_acc]=calc_xyar(Ne1(ii,jj),Ni1(ii,jj),Btot1(ii,jj));
    [eps]=calc_eps1(xar1,xar2,yar1,yar2);
    [nperps,nperpf,D]=calc_nperpsf(eps,npar);

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
% [dpsidr1]=calc_interp22(imin,imax,jmin,jmax,dpsidr);
% [dpsidz1]=calc_interp22(imin,imax,jmin,jmax,dpsidz);
gradpsi(ii,jj)=sqrt(dpsidr1(ii,jj).^2+dpsidz1(ii,jj).^2);
%cnr=(cnteta1(1,j)*dpsidz1(ii,jj)-cirho*cnrho*dpsidr1(ii,jj))/gradpsi(ii,jj);
%cnz=(-cnteta1(1,j)*dpsidr1(ii,jj)-cirho*cnrho*dpsidz1(ii,jj))/gradpsi(ii,jj);
cnz=(cnteta1(1,j)*Bz1(ii,jj)+cnrho*Br1(ii,jj))/sqrt(Br1(ii,jj)^2+Bz1(ii,jj)^2);
cnr=(cnteta1(1,j)*Br1(ii,jj)-cnrho*Bz1(ii,jj))/sqrt(Br1(ii,jj)^2+Bz1(ii,jj)^2);
cm=cnphi*eqdsk_r1(jj);

cnteta111=(cnr*Br1(ii,jj)+cnz*Bz1(ii,jj))/sqrt(Br1(ii,jj)^2+Bz1(ii,jj)^2);

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
ttstep(1,j)=0;      %error
pscatter1(1,j)=0;   %error
DD(1,j)=0;



%%
for i=1:lengthray
    
   
%     [ddd,ddt,ii,jj]=ddd1(z,r,phi,cnz,cnr,cm);
%     dddd(i,j)=ddd;
%     ddt1(i,j)=ddt;
%     ii1(i,j)=ii;
%     jj1(i,j)=jj;
%     sne1(i,j)=Ne1(ii,jj);
%     Btot(i,j)=sqrt(Br1(ii,jj)^2+Bz1(ii,jj)^2+Bt1(ii,jj)^2);
    
%% 
%if (Ne1(ii,jj)==0 || r<eqdsk_r(1) || r>eqdsk_r(end) || ddt.z>c || ddt.r>c) %ddt>c ï¿½ï¿½ï¿½ï¿½Ê²Ã´ï¿½ï¿½ï¿½â£¿
%if 1% inpolygon(r,z,gvar.rlim,gvar.zlim) &&
%inpolygon(r,z,gvar.rbbbs,gvar.zbbbs)==0 % just include SOL

if 1%rhopsi11(ii,jj)>0.9 && inpolygon(r,z,gvar.rlim,gvar.zlim)  % change the calcation zone
   % if inpolygon(r,z,gvar.rlim,gvar.zlim)
    [ddd,ddt,ii,jj]=calc_ddd1_old(z,r,phi,cnz,cnr,cm);
    %if abs(ddt.z)<c & abs(ddt.r)<c
    dddd2(i,j)=ddd;
    dddd1.cnz(i,j)=ddd.cnz;
    dddd1.cnr(i,j)=ddd.cnr;
    dddd1.cm(i,j)=ddd.cm;
    dddd1.w(i,j)=ddd.w;
    dddd1.r(i,j)=ddd.r;
    dddd1.z(i,j)=ddd.z;
    dddd1.phi(i,j)=ddd.phi;
   
    
    ddt1.cnz(i,j)=ddt.cnz;
    ddt1.cnr(i,j)=ddt.cnr;
    ddt1.cm(i,j)=ddt.cm;
    ddt1.r(i,j)=ddt.r;
    ddt1.z(i,j)=ddt.z;
    ddt1.phi(i,j)=ddt.phi;
    
    ii1(i,j)=ii;
    jj1(i,j)=jj;

    sne1(i,j)=Ne1(ii,jj);
    Btot(i,j)=sqrt(Br1(ii,jj)^2+Bz1(ii,jj)^2+Bt1(ii,jj)^2);
    [xar1,xar2,yar1,yar2,npar_acc]=calc_xyar(sne1(i,j),sne1(i,j),Btot(i,j));
    [eps]=calc_eps1(xar1,xar2,yar1,yar2);
    [nperps,nperpf,D]=calc_nperpsf(eps,Npar1(i,j));
    npar_acc1(i,j)=npar_acc;
    DD(i,j)=D;
    sbr(i,j)=Br1(ii,jj);
    sbz(i,j)=Bz1(ii,jj);
    sbt(i,j)=Bt1(ii,jj);
    nperps1(i,j)=imag(nperps);
    nperps11(i,j)=nperps;
    nperpf1(i,j)=imag(nperpf);
    nperpf11(i,j)=nperpf;
    ste(i,j)=Te1(ii,jj);
    sti(i,j)=Te1(ii,jj);
    
[cnprim]=absorb(nperps11(i,j),Npar1(i,j),ste(i,j),sti(i,j),sne1(i,j),sne1(i,j),Btot(i,j),Zeff,w,dddd2(i,j));
    cnprim1.e(i,j)=cnprim.e; 
    cnprim1.i(i,j)=cnprim.i;
    cnprim1.cl(i,j)=cnprim.cl;
%%
[z,r,phi,cnz,cnr,cm,us]=rk2_old(z,r,phi,cnz,cnr,cm,us);
%[ii,jj]=findposition(z,r);
[~,ii]=min(abs(z-eqdsk_z1));
[~,jj]=min(abs(r-eqdsk_r1));
% [imax,imin,jmax,jmin]=calc_position(z,r,eqdsk_z,eqdsk_r);
% calc_position_new;
z1(i+1,j)=z;
r1(i+1,j)=r;
phi1(i+1,j)=phi;
cnz1(i+1,j)=cnz;
cnr1(i+1,j)=cnr;
cm1(i+1,j)=cm;
us1(i+1,j)=us;
ttstep(i+1,j)=(z1(i+1,j)-z1(i,j))./ddt1.z(i,j);


[teta]=calc_teta11(z1(i+1,j),r1(i+1,j),z1(i,j),r1(i,j),gvar);
teta1(i+1,j)=teta1(i,j)+teta;
% if teta1(i,j)>teta
%     teta1(i+1,j)=teta+2*pi;
% else
%     teta1(i+1,j)=teta;
% end

Npar1(i+1,j)=(Br1(ii,jj)*cnr+Bz1(ii,jj)*cnz+Bt1(ii,jj)*cm/eqdsk_r1(jj))/Btot(i,j);
%


%%  ï¿½ï¿½ï¿½Ö¹ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½Ú²ï¿½ï¿½ï¿½ï¿½ï¿½Ë®ï¿½?
if  j~=1
    Deltan11(i,j)=Deltan1(ii,jj);
    if  Deltan1(ii,jj)~=0  &&  (imag(nperps)~=0 && imag(nperpf)~=0)~=1
        [pscatter,snperm,ds,cnz,cnr,cm,lent,lenl,beta1,Vgperm,nu1,ps]=densfluc_old(z,r,phi,cnz,cnr,cm,tstep,ktheta11,kmean);       %ï¿½ï¿½ï¿½ï¿½ï¿½ï²¢Ã»ï¿½Ð¸ï¿½ï¿½ï¿½ï¿½ÂµÄ£ï¿½cnz,cnr,cmï¿½ï¿½,ï¿½ï¿½ï¿½Þ¸ï¿½
         lent1(i,j)=lent;
         lenl1(i,j)=lenl;
         beta11(i,j)=beta1;
         Vgperm1(i,j)=Vgperm;
         nu11(i,j)=nu1;
        snperm1(i,j)=snperm;
        ps1(i,j)=ps;
        pscatter1(i,j)=pscatter;
       
    end
end

    if ((abs(Npar1(i,j))<npar_acc1(i,j))||(imag(nperps)>0&&imag(nperpf)>0))     %ï¿½ï¿½Ö¹ï¿½ï¿½ï¿½ï¿½
%         Npar1=Npar(i,j);
%         z1(i,j)=z;
%         r1(i,j)=r;
%         phi1(i,j)=phi;
%         us1(i,j)=us;c 
%         cnz1(i,j)=cnz;
%         cnr1(i,j)=cnr;
%         cm1(i,j)=cm;
%         rhopsi1(i,j)=rhopsi11(ii,jj);
%         cnteta(i,j)=(cnr*Br1(ii,jj)+cnz*Bz1(ii,jj))/sqrt(Br1(ii,jj)^2+Bz1(ii,jj)^2);
    %[cnprim]=absorb(nperps11(i,j),Npar1(i,j),ste(i,j),sti(i,j),sne1(i,j),sne1(i,j),Btot(i,j),Zeff,w);
%     cnprim1.e(i,j)=cnprim.e; 
%     cnprim1.i(i,j)=cnprim.i;
%     cnprim1.cl(i,j)=cnprim.cl;
%     deltaP.e(i,j)=exp(-2*cnprim1.e(i,j)*tstep);
%     deltaP.i(i,j)=exp(-2*cnprim1.i(i,j)*tstep);
%     deltaP.cl(i,j)=exp(-2*cnprim1.cl(i,j)*tstep);
        rhopsi1(i+1,j)=rhopsi11(ii,jj);
        cnteta1(i+1,j)=(cnr*Br1(ii,jj)+cnz*Bz1(ii,jj))/sqrt(Br1(ii,jj)^2+Bz1(ii,jj)^2);
%         break;
    else
        if rhopsi11(ii,jj)>1.2      %ï¿½ï¿½rhoÎª1.2ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½
            cnteta=(cnr*Br1(ii,jj)+cnz*Bz1(ii,jj))/sqrt(Br1(ii,jj)^2+Bz1(ii,jj)^2);
            cnpsi=(cnz*Br1(ii,jj)-cnr*Bz1(ii,jj))/sqrt(Br1(ii,jj)^2+Bz1(ii,jj)^2);
            cnzref=(cnteta*Bz1(ii,jj)-cnpsi*Br1(ii,jj))/sqrt(Br1(ii,jj)^2+Bz1(ii,jj)^2);
            cnrref=(cnteta*Br1(ii,jj)+cnpsi*Bz1(ii,jj))/sqrt(Br1(ii,jj)^2+Bz1(ii,jj)^2);
            cnz=cnzref;
            cnr=cnrref;
        end
        rhopsi1(i+1,j)=rhopsi11(ii,jj);
        cnteta1(i+1,j)=(cnr*Br1(ii,jj)+cnz*Bz1(ii,jj))/sqrt(Br1(ii,jj)^2+Bz1(ii,jj)^2);
    end
    end
end
%%
end

data_save;

toc
%figure
