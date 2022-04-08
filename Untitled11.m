global Bt1 Br1 Bz1 eqdsk_r1 eqdsk_z1  Ne1 Te1 w Btot1 Ni1 deltan1
    npar=-2.04;
    a=0.45;
    ksi0roui=0.1;
    cirho=1;

    z0=0;
    r0=2.33;
    phi0=0;
    [ii,jj]=findposition(z0,r0);
    Btot(ii,jj)=(Br1(ii,jj)^2+Bz1(ii,jj)^2+Bt1(ii,jj)^2)^0.5;
 [xar1,xar2,yar1,yar2,npar_acc]=calc_xyar(Ne1(ii,jj),Ni1(ii,jj),Btot1(ii,jj));
 [eps]=calc_eps1(xar1,xar2,yar1,yar2);
 [nperps,nperpf,D]=calc_nperpsf(eps,npar);
 %[eps,nperps,nperpf,D]=calc_eps(xar1,xar2,yar1,yar2,npar);
 kpar=w*npar./c;

% 
 if imag(nperps)==0
     nperp=nperps;
 end

  gradpsi=sqrt(dpsidr.^2+dpsidz.^2);
 
 cnteta1=0;
 cnphi=npar*Btot1(ii,jj)/Bt1(ii,jj);
 cn2=npar^2+nperps^2;
 cnrho2=cn2-cnteta1^2-cnphi^2;
 cnrho=sqrt(cnrho2);
 cnr0=(cirho*cnrho*dpsidr(ii,jj))/gradpsi(ii,jj);
 cnz0=(cirho*cnrho*dpsidz(ii,jj))/gradpsi(ii,jj);
 cm0=cnphi*eqdsk_r1(jj);
 
npar1=(cnr0*Br1(ii,jj)+cnz0*Bz1(ii,jj)+Bt1(ii,jj)*cnphi)/Btot1(ii,jj);
us=0;
 Btot2=[];
    nperps=[];
    nperpf=[];
    cnteta1=[];
    z1=[];
    r1=[];
    phi1=[];
    cnz1=[];
    cnr1=[];
    cm1=[];
    us1=[];
    Npar1=[];
    teta1=[];
     rhopsi1=[];    
    dddd1=[];
    ddd=[];
    ddt=[];
    ddt1=[];
    ii1=[];
    jj1=[];
    sne1=[];
    npar_acc1=[];
    sbr=[];
    sbz=[];
    sbt=[];
    nperps1=[];
    nperpf1=[];
    nperp=[]; 
    deltan11=[];
    xrandom1=[];
    snperm1=[];
    sbr1=[];
    sbz1=[];
    sbt1=[];
    nperps1=[];
    nperpf1=[];
    k=[];
    
parfor j=1:nray

k=[k,j];
%%

%  z1(1,j)=z;
%  r1(1,j)=r;
%  phi1(1,j)=phi;
%  cnz1(1,j)=cnz;
%  cnr1(1,j)=cnr;
%  cm1(1,j)=cm;
%  us1(1,j)=us;
%  Npar1(1,j)=npar1;     
%  teta1(1,j)=0;
    z=z0;
    r=r0;
    phi=phi0;
    cnz=cnz0;
    cnr=cnr0;
    cm=cm0;
    us=us0;
    [a]=solve1(r);
    [r1]=[r1,a];
end
%[ddd,ddt,ii,jj]=calc_ddd1(z,r,phi,cnz,cnr,cm);
 %for i=1:lengthray
    % complex(i,j)

%  if inpolygon(r,z,gvar.rlim,gvar.zlim) & inpolygon(r,z,gvar.rbbbs,gvar.zbbbs)==0
% % %if inpolygon(r,z,gvar.rlim,gvar.zlim)


%        dddd1=[dddd1;ddd];
%       ddt1=[ddt1;ddt];
%       ii1=[ii1;ii];
%       jj1=[jj1;jj];
% %      if abs(ddt.z)<c && abs(ddt.r)<c
% %     end
% 
%       sne=Ne1(ii,jj);
%       Btot=sqrt(Br1(ii,jj)^2+Bz1(ii,jj)^2+Bt1(ii,jj)^2);
%      [xar1,xar2,yar1,yar2,npar_acc]=calc_xyar(sne,sne,Btot);
%      [eps]=calc_eps1(xar1,xar2,yar1,yar2);
%      [nperps,nperpf,D]=calc_nperpsf(eps,npar1);
%     
%      npar_acc1=[npar_acc1,npar_acc];
%      sbr=Br1(ii,jj);
%      sbz=Bz1(ii,jj);
%      sbt=Bt1(ii,jj);
%      sbr1=[sbr1;sbr];
%      sbz1=[sbz1;sbr];
%      sbt1=[sbt1;sbr];
%      nperps2=imag(nperps);
%      nperpf2=imag(nperpf);
%      nperps1=[nperps1;nperps2];
%      nperpf1=[nperpf1;nperpf2];
% 
  % [z,r,phi,cnz,cnr,cm,us]=rk2(z,r,phi,cnz,cnr,cm,us);
% %  [ii,jj]=findposition(z,r);
%  z1=[z1;z]
%  
%    r1=[r1;r];
%   phi1=[phi1;phi];
%  cnz1=[cnz1;cnz];
%  cnr1=[cnr1;cnr];
%  cm1=[cm1;cm];
%  us1=[us1;us];

% [teta]=calc_teta11(z1(i+1,j),r1(i+1,j),z1(i,j),r1(i,j),gvar);
% teta1(i+1,j)=teta1(i,j)+teta;
% % if teta1(i,j)>teta
% %     teta1(i+1,j)=teta+2*pi;
% % else
% %     teta1(i+1,j)=teta;

% end
%  z1=[z1,z]
%  r1=[r1,r];
%  
% end
% Npar1(i+1,j)=(Br1(ii,jj)*cnr+Bz1(ii,jj)*cnz+Bt1(ii,jj)*cm/eqdsk_r1(jj))/Btot(i,j);
% % z1(i+1,j)=z+tstep*ddt.z;
% % r1(i+1,j)=r+tstep*ddt.r;
% % phi1(i+1,j)=phi+tstep*ddt.phi;
% % cnz1(i+1,j)=cnz+tstep*ddt.cnz;
% % cnr1(i+1,j)=cnr+tstep*ddt.cnr;
% % cm1(i+1,j)=cm+tstep*ddt.cm;
% % us1(i+1,j)=us+sqrt((tstep*ddt.z)^2+(tstep*ddt.r)^2+(tstep*ddt.phi)^2);
% % z=z1(i+1,j);
% % r=r1(i+1,j);
% % phi=phi1(i+1,j);
% % cnz=cnz1(i+1,j);
% % cnr=cnr1(i+1,j);
% % cm=cm1(i+1,j);
% % us=us1(i+1,j);
% 
% %%  区分刮削层和内部涨落水平
% if  j~=1
%     deltan11(i,j)=deltan1(ii,jj);
%     if  deltan1(ii,jj)~=0  &&  (imag(nperps)~=0 && imag(nperpf)~=0)~=1
%         [snperm,ds,cnz,cnr,cm,xrandom]=densfluc(z,r,phi,cnz,cnr,cm,tstep);       %在这里并没有给出新的（cnz,cnr,cm）,已修改
%         xrandom1(i,j)=xrandom;
%         snperm1(i,j)=snperm;
%        
%     end
% % else
% %     ds=0;
% end
% %  us=us+ds;      %意义在哪？如果(z,r)不变的话，那么ds的变化是怎么来的？
% 
% %% 
% %[z,r,phi,cnz,cnr,cm,us]=rk11(z,r,phi,cnz,cnr,cm,us);
% %[z,r,phi,cnz,cnr,cm,us]=rk2(z,r,phi,cnz,cnr,cm,us);
% 
% 
% 
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
%     end
% end
% %%
% end
% end