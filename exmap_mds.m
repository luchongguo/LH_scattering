function [npsip,nrhot,rhotb,Rin,Rout,rg]=exmap_mds(R,Z,t,shot,icheck)
%********************************************************
%  mapping to flux coordinate
%
%    input:
%       (R,Z), postions of the measurments,
%           t, the time slices
%         shot,  shot number,
%     if icheck=1, plot the mapping of the surface
%                 at each time slice and pause for check.
%
%    output:
%       npsip, normlaized ploidal flux,
%      nrhot, normalized rho_t=sqrt(psit/B0/pi), 
%              psit is the toroidal flux,
%      rhotb,  rhotb=rhot(npsit=1), i.e. at boundary
%   [Rin,Rout], the mapping major radius
%               at high and low field side at Z=Zc 
%
%     by Y. Sun
%**************************************************** ****
%addpath('../mds/');
%----- load eq data
   [t0,rg,zg,npsi,Rc,Zc,gpsip,nrhot0,rhotb0]=geteq(shot);
 
%-----
it0=0;
   for i=1:length(t)
       [dum,it]=min(abs(t0-t(i)));
        tit=['@t=',num2str(t0(it)),'s'];
        if it~=it0
        [npsip(:,i),nrhot(:,i),Rin(:,i),Rout(:,i)]= ... 
           calrho(rg,zg,npsi(:,:,it),Rc(it),Zc(it),gpsip,nrhot0(:,it),R,Z,tit,icheck);
        else
           npsip(:,i)=npsip(:,i-1);
           nrhot(:,i)=nrhot(:,i-1);
           Rin(:,i)=Rin(:,i-1);   
           Rout(:,i)=Rout(:,i-1);
        end
       rhotb(i)=rhotb0(it);
       it0=it;
   end           
end

function [npsip,nrhot,Rin,Rout]=calrho(rg,zg,npsi,Rc,Zc,gpsip,nrhot,R,Z,tit,icheck)
%***************************************************************
%   cal rho
%  INPUTS
%      npsi(rg,zg), the 2D normalized poloidal flux profile
%                   on the grids [rg, zg]
%      [Rc,Zc], position of the magnetic axis
%      gpsip, equaly spaced normalized poloidal flux grids
%      nrhot(gpsip), the normalized sqrt(\psi_t)
%      (R,Z), postion of the measurments,
%       tit, the time info in the title in the checking figure
%       icheck, if =1, plot the positions 
%                    superimposed by the flux surfaces   
%  OUTputs:   
%       npsip, normlaized ploidal flux,
%      nrhot, normalized rho_t=sqrt(psit/B0/pi), 
%              psit is the toroidal flux,
%      rhotb,  rhotb=rhot(npsit=1), i.e. at boundary
%   [Rin,Rout], the mapping major radius
%               at high and low field side at Z=Zc        
%***************************************************************
 %----- mapping psi_p  
   npsip=interp2(rg,zg,npsi,R,Z,'cubic');
%----- mapping Rout
   [Rin,Rout]=calR(rg,zg,npsi,Rc,Zc,npsip);
%----- mapping psi_t, rho_t         
  nrhot=interp1(sqrt(gpsip),nrhot,sqrt(npsip),'linear','extrap');    
%-----  icheck mapping
         if icheck==1
            myplot(rg,zg,npsi,R,Z,tit);
           disp(['   R(m)',',     Z(m)',',    psi_p',',     R_in',',    R_out'])
           [R',Z',npsip',Rin',Rout']
           pause
         end
end


function [t0,rg,zg,npsi,Rc,Zc,gpsip,nrhot,rhotb]=geteq(shot)
%***************************************************************
%   get the equilibrium profile from the gfile
%      npsi(rg,zg), the 2D normalized poloidal flux profile
%                   on the grids [rg, zg]
%      [Rc,Zc], position of the magnetic axis
%      gpsip, equaly spaced normalized poloidal flux grids
%      nrhot(gpsip), the normalized sqrt(\psi_t)
%      rhotb = rhot(psi=psi_b)
%      t0, the t slices
%      q(gpsip), the q profile
%      B0, toroidal magnetic field at the geometric center
%***************************************************************
 mdsconnect('202.127.204.12')
 shotnumber=mdsopen('efit_east',shot)
   rg=mdsvalue2(['\r']);
   zg=mdsvalue2(['\z']);
  B0=mdsvalue2(['\BCENTR'])'; 
  t0=mdsvalue2(['dim_of(\q95)']);
    nt= length(t0);
%----- q 
      q=mdsvalue2(['\qpsi']);  
      nw=length(q)/nt;
  gpsip=linspace(0,1,nw);
  q=reshape(q,[nw,nt]);
 %---- psi_p  psit rho_t   
    psimag=mdsvalue2(['\SSIMAG']);  
    psibry=mdsvalue2(['\SSIBRY']); 
     psibm=psibry-psimag;
  psirz=mdsvalue2(['\psirz']);
  psirz=reshape(psirz,[length(rg) length(zg) nt]);
    for i=1:nt
        psirz(:,:,i)=  psirz(:,:,i)';
        npsi(:,:,i)=(psirz(:,:,i)-psimag(i))./(psibry(i)-psimag(i));
           psip0=psibm(i)*gpsip;                 % 1-D psi_p in real unit.
           psit(:,i)=calpsit(psip0,q(:,i));      % 1-D psi_t in real unit.
          rhot(:,i)=sqrt(abs(2*psit(:,i)/B0(i))); % 1-D rho_t in real unit. Note the psi in MDS is normalized to 2*pi.
        rhotb(i)=rhot(end,i);
        nrhot(:,i)=rhot(:,i)/rhotb(i);
    end
 %---- psi_p  
    Rc=mdsvalue2(['\RMAXIS']);  
    Zc=mdsvalue2(['\ZMAXIS']); 
end

function a=mdsvalue2(sig)
%**********************************************************************
%   use reshape to change it into 1 collumn
%   this is necessary for the new data type saved in the server.
%**********************************************************************
a=mdsvalue(sig); 
[m,n]=size(a);
a=reshape(a,m*n,1);
end

function myplot(rg,zg,npsi,R,Z,tit)
%***************************************************************
%  plot the postions for mapping with plasma configuration
%***************************************************************
           figure(111)
           clf
           contour(rg,zg,npsi,20);
           hold on
          contour(rg,zg,npsi,[1 1],'--r');  
           plot(R,Z,'o','linewidth',1.5);
            xlabel('R(m)','fontsize',15)
            ylabel('Z(m)','fontsize',15)
            title(['Location of the measurment',tit],'fontsize',15)
                set(gca,'fontsize',12)
           hold off
           axis image
end

function [Rin,Rout]=calR(rg,zg,psirz,Rc,Zc,psip)
%*********************************
% calculate R_in and R_out
%    at the midplane
%*********************************
nR=101;
Rin0=linspace(rg(1),Rc,nR);
Rout0=linspace(Rc,rg(end),nR);
Z0=Zc*ones(1,nR);
psi_in=interp2(rg,zg,psirz,Rin0,Z0,'cubic');
psi_out=interp2(rg,zg,psirz,Rout0,Z0,'cubic');
Rin=interp1(psi_in,Rin0,psip,'cubic');
Rout=interp1(psi_out,Rout0,psip,'cubic');
end

function psit=calpsit(psip,q)
%*********************************
% calculate psit
%*********************************
lp=length(psip);
psipm=mean([psip(1:(lp-1)); psip(2:lp)]);
nq=interp1(psip,q,psipm,'cubic');
dpsip=diff(psip);
psit=cumsum(nq.*dpsip);
psit=[0 psit];
end