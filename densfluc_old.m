function [pscatter,snperm,ds,cnz,cnr,cm,lent,lenl,beta1,Vgperm,nu1,ps]=densfluc_old(z,r,phi,cnz,cnr,cm,tstep,ktheta11,kmean)
%function [cnz,cnr,cm]=densfluc()
% density fluction
%input: (z,r,phi),(cnz,cnr,cm),B,gamma;
%output: (cnz,cnr,cm)
%%
global Bt1 Br1 Bz1 eqdsk_r1 eqdsk_z1 Ne1 Ni1 w Deltan1 Ktheta1
%global Bt Bp Br Bz eqdsk_r eqdsk_z  Ne Ni Te w Btot Deltan rhopsi
% z=ray.z;
% r=ray.r;
% phi=ray.phi;
% cnz=ray.cnz;
% cnr=ray.cnr;
% cm=ray.cm;
% tstep=1e-11;
[~,ii]=min(abs(z-eqdsk_z1));
[~,jj]=min(abs(r-eqdsk_r1));
% [imax,imin,jmax,jmin]=calc_position(z,r,eqdsk_z,eqdsk_r);
% calc_position_new;

br=Br1(ii,jj);
bz=Bz1(ii,jj);
bphi=Bt1(ii,jj);
ne=Ne1(ii,jj);
ni=Ni1(ii,jj);
deltan=Deltan1(ii,jj);
%ktheta=abs(2*pi*100*Ktheta1(ii,jj));  %��ôȷ����
%ktheta=abs(2*pi*100);
if ktheta11==0
    ktheta11=1e-20;
end
ktheta=200*pi/abs(ktheta11);
pscatter=0;
c=299792458;

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
    [lent,lenl,len,cnzfs,cnrfs,cmfs,beta1,Vgperm,nu1]=Copy_of_ks(z,r,phi,cnz,cnr,cm,cnpar,skperp,skperm,ktheta,deltan,kmean);
    lenfs=lenl;                      %  scattering length
    if(ratk<1)
        pfs=lenfs*ktheta;
    else
        pfs=lenfs*ktheta/ratk;
    end
    ps=pfs;
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
     [lent,lenl,len,cnzff,cnrff,cmff,beta1,Vgperm,nu1]=Copy_of_ks(z,r,phi,cnz,cnr,cm,cnpar,skperp,skperm,ktheta,deltan,kmean);
    lenff=lenl;                     % scattering length
    if(ratk<1)
        pff=lenff*ktheta;
    else
        pff=lenff*ktheta/ratk;
    end
    ps=pff;
    if(pff<1)
        ds=0;
        xrandom=0;
        return;
    end
    sumff=lent;
    
    %deltt=tavg/(sumff+sumfs);
    deltt=tstep;
    [~,ddt]=calc_ddd1_old(z,r,phi,cnz,cnr,cm);     %dddz,dddr,dddphi��,��ʱ������ɢ����ֵ��
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
    probfs=deltt*sumfs;        % FW->SW ����        �������⣬deltt�Ǿ������ʱ�䣬������Ӧ����ʱ������
    probff=deltt*sumff;        % FW->FW ����
    xrandom=rand(1);
    probf=probff+probfs;
    
    pscatter=probf;  %20210311
    
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
     [lent,lenl,len,cnzss,cnrss,cmss,beta1,Vgperm,nu1]=Copy_of_ks(z,r,phi,cnz,cnr,cm,cnpar,skperp,skperm,ktheta,deltan,kmean);
    lenfs=lenl;                      %  scattering length
    if(ratk<1)
        pfs=lenfs*ktheta;
    else
        pfs=lenfs*ktheta/ratk;
    end
    ps=pfs;
    if(pfs<1 || pfs==1) || isnan(pfs)
        ds=0;
        xrandom=0;
        return;
    end
    sumss=lent;
    
    %% scattering condition
    %deltt=tavg/sumss;
    deltt=tstep;
    [~,ddt]=calc_ddd1_old(z,r,phi,cnz,cnr,cm);      %dddz,dddr,dddphi��,��ʱ�Ѿ���ɢ����ֵ��
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
    %probss=deltt*sumss;        % SW->SW ����
    probss=deltt*sumss;
    pscatter=probss;
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
     [lent,lenl,len,cnzsf,cnrsf,cmsf,beta1,Vgperm,nu1]=Copy_of_ks(z,r,phi,cnz,cnr,cm,cnpar,skperp,skperm,ktheta,deltan,kmean);
    lensf=lenl;                      %  scattering length
    if(ratk<1)
        psf=lensf*ktheta;
    else
        psf=lensf*ktheta/ratk;
    end
    ps=psf;
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
    [lent,lenl,len,cnzss,cnrss,cmss,beta1,Vgperm,nu1]=Copy_of_ks(z,r,phi,cnz,cnr,cm,cnpar,skperp,skperm,ktheta,deltan,kmean);
    lenss=lenl;                     % scattering length
    if(ratk<1)
        pss=lenss*ktheta;
    else
        pss=lenss*ktheta/ratk;
    end
    ps=pss;
    if(pss<1)
        ds=0;
        xrandom=0;
        return;
    end
    sumss=lent;
    
    %deltt=tavg/(sumss+sumsf);
    deltt=tstep;
    [~,ddt]=calc_ddd1_old(z,r,phi,cnz,cnr,cm);     %dddz,dddr,dddphi��,��ʱ�Ѿ���ɢ����ֵ��
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
    probsf=deltt*sumsf;        % SW->FW ����
    probss=deltt*sumss;        % SW->SW ����
    xrandom=rand(1);
    probs=probss+probsf;
    pscatter=probs;
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