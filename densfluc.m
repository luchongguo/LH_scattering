function [snperm,ds,cnz,cnr,cm,xrandom]=densfluc(z,r,phi,cnz,cnr,cm,tstep)
%function [cnz,cnr,cm]=densfluc()
% density fluction
%input: (z,r,phi),(cnz,cnr,cm),B,gamma;
%output: (cnz,cnr,cm)
%%
%global Bt1 Br1 Bz1 eqdsk_r1 eqdsk_z1 Ne1 Ni1 w Deltan1
global Bt Bp Br Bz eqdsk_r eqdsk_z  Ne Ni Te w Btot Deltan rhopsi
% z=ray.z;
% r=ray.r;
% phi=ray.phi;
% cnz=ray.cnz;
% cnr=ray.cnr;
% cm=ray.cm;
% tstep=1e-11;
% [~,ii]=min(abs(z-eqdsk_z1));
% [~,jj]=min(abs(r-eqdsk_r1));
[imax,imin,jmax,jmin]=calc_position(z,r,eqdsk_z,eqdsk_r);
calc_position_new;

br=Br1(ii,jj);
bz=Bz1(ii,jj);
bphi=Bt1(ii,jj);
ne=Ne1(ii,jj);
ni=Ni1(ii,jj);
deltan=Deltan1(ii,jj);

% dne=dNe1(ii,jj);
% dbrr=dBrr1(ii,jj);
% dbzz=dBzz1(ii,jj);

c=299792458;
ktheta=2*pi*100;        %怎么确定？
tavg=1e-2;
%%      initial mode of LHW (FW | SW)
[gamma,cn]=calc_gamma(bz,br,bphi,cnz,cnr,cm,r);
btot=sqrt(bphi^2+bz^2+br^2);
cnperp=cn*sin(gamma);
cnpar=cn*cos(gamma);
[xar1,xar2,yar1,yar2]=calc_xyar(ne,ni,btot);
[eps]=calc_eps1(xar1,xar2,yar1,yar2);
[nperps,nperpf,~]=calc_nperpsf(eps,cnpar);
%[~,nperps,nperpf]=calc_eps(xar1,xar2,yar1,yar2,cnpar);
cnperps=nperps;
cnperpf=nperpf;

difcnperp1=abs(cnperp-cnperps);
difcnperp2=abs(cnperp-cnperpf);
if difcnperp1>difcnperp2        %compare the nperp in order to find the mode of LHW
    mode=1;                     %FW
else
    mode=2;                     %SW
end
%%  FW->SW   becasue of fluctuation
if mode==1
    
    snperp=cnperpf;             % Nperp before scattering 
    snperm=cnperps;             % Nperp after scattering
    skperp=snperp*w/c;          % kperp before sacttering
    skperm=snperm*w/c;
    
    ratk=(skperp/ktheta)^2;      % (kperp/ktheta)^2
    [lent,lenl,len,cnzfs,cnrfs,cmfs]=ks(z,r,phi,cnz,cnr,cm,cnpar,skperp,skperm,ktheta,deltan);
    lenfs=lenl;                      %  scattering length
    if(ratk<1)
        pfs=lenfs*ktheta;
    else
        pfs=lenfs*ktheta/ratk;
    end
    if(pfs<1&pfs==1) || isnan(pfs)
        ds=0;
        xrandom=0;
        return;
    end
    sumfs=lent;
    %%  FW->FW   becasue of fluctuation
    snperp=cnperpf;             % Nperp before scattering 
    snperm=cnperpf;             % Nperp after scattering
    skperp=snperp*w/c;          % kperp before sacttering
    skperm=snperm*w/c; 
   
    ratk=(skperp/ktheta)^2;      % (kperp/ktheta)^2
     [lent,lenl,len,cnzff,cnrff,cmff]=ks(z,r,phi,cnz,cnr,cm,cnpar,skperp,skperm,ktheta,deltan);
    lenff=lenl;                     % scattering length
    if(ratk<1)
        pff=lenff*ktheta;
    else
        pff=lenff*ktheta/ratk;
    end
    if(pff<1)
        ds=0;
        xrandom=0;
        return;
    end
    sumff=lent;
    
    %deltt=tavg/(sumff+sumfs);
    deltt=tstep;
    [~,ddt]=calc_ddd1(z,r,phi,cnz,cnr,cm);     %dddz,dddr,dddphi等,此时还不是散射后的值了
    vg=(ddt.z^2+ddt.r^2+(r*ddt.phi)^2)^0.5;
    ds=vg*deltt;
%     sout=sin1+ds;
%     deltast=stout-stin;
%     deltat=min(deltast,ds);
%     if (stout<sout)
%         return;
%     elseif (stout>sout)
%         sin1=sout;
%         stout=stin+deltast;
%         stin=sin1;
%     end
    probfs=deltt*sumfs;        % FW->SW 概率        存在问题，deltt是距离而非时间，但书上应该是时间间隔？
    probff=deltt*sumff;        % FW->FW 概率
    xrandom=rand(1);
    probf=probff+probfs;
    if (xrandom>probf)                % no scattering
        snperm=cnperpf;
    elseif (xrandom>probff & xrandom<probf) % FW->SW
        cnz=cnzfs;
        cnr=cnrfs;
        cm=cmfs;
        snperm=cnperps;  
    elseif (xrandom<probff)
        cnz=cnzff;
        cnr=cnrff;
        cm=cmff;
        snperm=cnperpf;
    end
end    
    

    
    
%%
%%  SW->SW   becasue of fluctuation
if mode==2 & imag(cnperpf)~=0  
    
    snperp=cnperps;             % Nperp before scattering 
    snperm=cnperps;             % Nperp after scattering
    skperp=snperp*w/c;          % kperp before sacttering
    skperm=snperp*w/c;
    ratk=(skperp/ktheta)^2;      % (kperp/ktheta)^2
     [lent,lenl,len,cnzss,cnrss,cmss]=ks(z,r,phi,cnz,cnr,cm,cnpar,skperp,skperm,ktheta,deltan);
    lenfs=lenl;                      %  scattering length
    if(ratk<1)
        pfs=lenfs*ktheta;
    else
        pfs=lenfs*ktheta/ratk;
    end
    if(pfs<1 || pfs==1) || isnan(pfs)
        ds=0;
        xrandom=0;
        return;
    end
    sumss=lent;
    
    %% scattering condition
    %deltt=tavg/sumss;
    deltt=tstep;
    [~,ddt]=calc_ddd1(z,r,phi,cnz,cnr,cm);      %dddz,dddr,dddphi等,此时已经是散射后的值了
    vg=sqrt(ddt.z^2+ddt.r^2+(r*ddt.phi)^2);
    ds=vg*deltt;
%     sin1=0;
%     sout=sin1+ds;
%     deltast=stout-stin;
%     deltat=min(deltast,ds);
%     if (stout<sout)
%         return;
%     elseif (stout>sout)
%         sin1=sout;
%         stout=stin+deltast;
%         stin=sin1;
%     end
    %probss=deltt*sumss;        % SW->SW 概率
    probss=deltt*sumss;
    xrandom=rand(1);
    if (xrandom>probss)               % no scattering
        
        return;
    else                        % SW->SW
        cnz=cnzss;
        cnr=cnrss;
        cm=cmss;
        snperm=cnperps; 
    end
end
    
if mode==2 & imag(cnperpf)==0
    %%  SW-> FW   becasue of fluctuation
    snperp=cnperps;             % Nperp before scattering 
    snperm=cnperpf;             % Nperp after scattering
    skperp=snperp*w/c;          % kperp before sacttering
    skperm=snperm*w/c;          % kperp after sacttering
    ratk=(skperp/ktheta)^2;      % (kperp/ktheta)^2
     [lent,lenl,len,cnzsf,cnrsf,cmsf]=ks(z,r,phi,cnz,cnr,cm,cnpar,skperp,skperm,ktheta,deltan);
    lensf=lenl;                      %  scattering length
    if(ratk<1)
        psf=lensf*ktheta;
    else
        psf=lensf*ktheta/ratk;
    end
    if(psf<1 || psf==1) || isnan(psf)
        ds=0;
        xrandom=0;
        return;
    end
    sumsf=lent;
    %%  SW->SW   becasue of fluctuation
    snperp=cnperps;             % Nperp before scattering 
    snperm=cnperps;             % Nperp after scattering
    skperp=snperp*w/c;          % kperp before sacttering
    ratk=(skperp/ktheta)^2;      % (kperp/ktheta)^2
    [lent,lenl,len,cnzss,cnrss,cmss]=ks(z,r,phi,cnz,cnr,cm,cnpar,skperp,skperm,ktheta,deltan);
    lenss=lenl;                     % scattering length
    if(ratk<1)
        pss=lenss*ktheta;
    else
        pss=lenss*ktheta/ratk;
    end
    if(pss<1)
        ds=0;
        xrandom=0;
        return;
    end
    sumss=lent;
    
    %deltt=tavg/(sumss+sumsf);
    deltt=tstep;
    [~,ddt]=calc_ddd1(z,r,phi,cnz,cnr,cm);     %dddz,dddr,dddphi等,此时已经是散射后的值了
    vg=(ddt.z^2+ddt.r^2+(r*ddt.phi)^2)^0.5;
    ds=vg*deltt;
%     sout=sin1+ds;
%     deltast=stout-stin;
%     deltat=min(deltast,ds);
%     if (stout<sout)
%         return;
%     elseif (stout>sout)
%         sin1=sout;
%         stout=stin+deltast;
%         stin=sin1;
%     end
    probsf=deltt*sumsf;        % SW->FW 概率
    probss=deltt*sumss;        % SW->SW 概率
    xrandom=rand(1);
    probs=probss+probsf;
    if (xrandom>probs)                % no scattering
        snperm=cnperps;
    elseif (xrandom>probss & xrandom<probs) % SW->FW
        cnz=cnzsf;
        cnr=cnrsf;
        cm=cmsf;
        snperm=cnperpf;    
    elseif (xrandom<probss)
        cnz=cnzss;
        cnr=cnrss;
        cm=cmss;
        snperm=cnperps;
    end 
end