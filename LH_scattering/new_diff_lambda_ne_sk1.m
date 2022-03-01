
clear
%clc
%npar1=[-1.82 -1.93 -2.04 -2.26 -2.48]; %4p6
%npar1=[-1.82  -2.23  -2.84]; %4p6
%npar1=[-1.82 -1.93 -2.1 -2.23 -2.37 -2.6 -2.84] %2p45
%npar1=[-1.82  -2.23 -2.84]
npar1=-2.23;
sktheta=2*pi*100/1;
kmean=0;
i_npar=length(npar1);

%w=4.6*2*pi*10^9; %低杂波波源频率

    ksi0roui=0.1;
    cirho=-1;
zz=0;rr1=[2.345:-0.01:2.295];

%zz=-0.3;rr1=[2.27296875000000,2.26203125000000,2.25000000000000,2.23687500000000]
%zz=0.3;rr1=[2.28390625000000,2.27296875000000,2.26203125000000,2.24890625000000];
phi=0;
%r=rr1(2)
r=2.34
z=zz;

%%
%load('F:\cbwu\办公\density_fluction_0408\testbeta\test_beta.mat')
%load('F:\cbwu\办公\density_fluction_0408\testbeta\data_original_f_2p45_90331_1.mat');
%load('F:\cbwu\办公\density_fluction_0408\testbeta\data_original_f_2p45_90331_5e18_lowdelta_0p5.mat');
%load('F:\cbwu\办公\density_fluction_0408\densityflucdata\90331_shenma\origin_data\data_f_2p45_90331_10e18_delta_0p2.mat');
load('E:\OneDrive - mail.ustc.edu.cn\个人信息\学籍相关\大论文\大论文数据\第五章\不同密度衰减长度\data_f_2p45_90331_8e18_lam_0p1_delta_0p2_new.mat');
%load('F:\cbwu\办公\density_fluction_0408\testbeta\data_original_f_2p45_90331_5e18_highdelta_0p1.mat');
color_ind=2;
tstep=0.5e-11;
[~,ii]=min(abs(z-eqdsk_z1));
[~,jj]=min(abs(r-eqdsk_r1));
br=Br1(ii,jj);
bz=Bz1(ii,jj);
bphi=Bt1(ii,jj);
ne=Ne1(ii,jj)
ni=Ne1(ii,jj);
deltan=Deltan1(ii,jj)
rhopsi11(ii,jj)
figure (1)
hold on
plot(rhopsi11(641,641:end),Ne(641,641:end),'-','linewidth',2)
plot(rhopsi11(641,347:end),Ne(641,347:end),'-','linewidth',2)

wpe=(ne*e*e/ee/m).^0.5;   %电子等离子体频率
wce=e*Btot1(ii,jj)/m;  %电子回旋频率
wpi=(ni*e*e/ee/M).^0.5;
wci=e*Btot1(ii,jj)/M;
npar_acc=wpe./wce+(1+(wpe./wce).^2-(wpi./w).^2).^(0.5); %可近性条件

xar1=(wpe./w).^2;
yar1=wce./w;
xar2=(wpi./w).^2;
yar2=wci./w;
% eps.perp=1+xar1/(yar1)^2-xar2;
% eps.par=1-xar1-xar2;
% eps.xy=xar1/yar1;
[eps]=calc_eps1(xar1,xar2,yar1,yar2);
%%
for i=1:1:i_npar
npar=npar1(i);
[nperps,nperpf,D]=calc_nperpsf(eps,npar);
j=1;
cnteta1(1,j)=0;
cnphi=npar*Btot1(ii,jj)/Bt1(ii,jj);
cn2=npar^2+nperps^2;
cnrho2=cn2-cnteta1(1,j)^2-cnphi^2;
cnrho=sqrt(cnrho2);
gradpsi(ii,jj)=sqrt(dpsidr1(ii,jj).^2+dpsidz1(ii,jj).^2);
% cnr=(cnteta1(1,j)*dpsidz1(ii,jj)-cirho*cnrho*dpsidr1(ii,jj))/gradpsi(ii,jj);
% cnz=(-cnteta1(1,j)*dpsidr1(ii,jj)-cirho*cnrho*dpsidz1(ii,jj))/gradpsi(ii,jj);
cnz=(cnteta1(1,j)*Bz1(ii,jj)+cnrho*Br1(ii,jj))/sqrt(Br1(ii,jj)^2+Bz1(ii,jj)^2);
cnr=(cnteta1(1,j)*Br1(ii,jj)-cnrho*Bz1(ii,jj))/sqrt(Br1(ii,jj)^2+Bz1(ii,jj)^2);
cm=cnphi*eqdsk_r1(jj);

uxx=eps.perp-1;
uxy=eps.xy;
uzz=eps.par-1;


%%      求H的三个分量，用于求体积V
hxx=1+wpe.^2./wce.^2;
hxy=0.5*wpe.^2./(w*wce);
hzz=1;
%%
skpr=npar*w/c;
skpr2=skpr.^2;

skperp=nperps*w/c;
skperm=nperps*w/c;
skperm1(i)=skperm;
skperp2=skperp.^2;
skperm2=skperm.^2;

dtor=(skpr2-eps.perp)*(skpr2+skperp2-eps.perp)-eps.xy^2;
dtorm=(skpr2-eps.perp)*(skpr2+skperm2-eps.perp)-eps.xy^2;

aa=skpr*skperp*(skpr2+skperp2-eps.perp)/dtor;
aam=skpr*skperm*(skpr2+skperm2-eps.perp)/dtorm;
bb=skpr*skperp*eps.xy/dtor;
bbm=skpr*skperm*eps.xy/dtorm;
aa1(i)=aa;
bb1(i)=bb;

u1=aa*aam+bb*bbm;
u2=aa*bbm+aam*bb;
urr=uxx*u1+uxy*u2;
uii=uxy*u1+uxx*u2;

Hh=hxx.*(aa.^2+bb.^2)+2.*aa.*bb.*hxy+hzz;
Hh1(i)=Hh;
Hhm=hxx.*(aam.^2+bbm.^2)+2.*aam.*bbm.*hxy+hzz;
ktot2=-skperp2*((wpe./wce).^2-(wpi./w).^2)-skpr2.*(-wpe.^2-wpi.^2)./w./w;
skpr21(i)=skpr2;
ktot21(i)=ktot2;
deltaw(i)=2./w*(skperp2./(ktot2)*(1+wpe.^2./wce.^2)+skpr2./(ktot2));
aaa(i)=2*pi*skperm./(deltaw(i).^2);

%%

beta3=-pi:0.01:pi;
% cosbeta1=(skperp^2+skperp^2-sktheta^2)./(2*skperp*skperp);
% beta2=acos(cosbeta1);
%  
for iii=1:length(beta3)
beta1(iii)=beta3(iii);

[cnzm,cnrm,cmm]=rotate(cnz,cnr,cm,r,br,bz,bphi,npar,skperp,skperm,w,beta1(iii)); % 计算新的（cnz,cnr,cm）下的群速度
%ddd %计算各个量对时间的求导
cnzm1(i,iii)=cnzm;
cnrm1(i,iii)=cnrm;
cmm1(i,iii)=cmm;
[~,ddt]=calc_ddd1_old(z,r,phi,cnzm,cnrm,cmm);
ddt1(i,iii)=ddt;
vgz=ddt.z;      %计算群速度的三个分量
vgr=ddt.r;
vgphi=r*ddt.phi;

% Vgdotb=vgz*bz+vgr*br+vgphi*bphi;
 btot2=bz^2+br^2+bphi^2;
% Vgperms=(vgz-Vgdotb*bz/btot)^2+(vgr-Vgdotb*br/btot)^2+(vgphi-Vgdotb*bphi/btot)^2
Vgtot2=vgz^2+vgr^2+vgphi^2;
Vgtot22(i,iii)=Vgtot2;
btot=sqrt(btot2);
Vgpar=(vgz*bz+vgr*br+vgphi*bphi)/btot;

Vgperms=(Vgtot2-Vgpar^2);       
Vgperm=sqrt(Vgperms);
Vgperm1(i,iii)=Vgperm;
cj11=uzz+cos(beta1(iii)).*urr;
cj12=sin(beta1(iii)).*uii;
cj1=cj11.^2+cj12.^2;



 lenls2(i,iii)=skperp^2+skperm^2-2*skperp*skperm*cos(beta1(iii));
 if beta1(iii)<0
     lenls(i,iii)=-sqrt(lenls2(i,iii));
 else
    lenls(i,iii)=sqrt(lenls2(i,iii));
 end
% if abs(lenls(i,iii))<abs(200*pi*kmean)
%     expbeta(i,iii)=0;
% else
expbeta(i,iii)=exp(-(lenls(i,iii)-200*pi*kmean).^2./(sktheta)^2);
%end
% P_ott1=(1+2*skperp2./ktot2.*((wpe./wce).^2-(wpi./w).^2).*(sin(beta1/2)).^2+(sin(beta1)).^2*wpe.^4/w/w/wce/wce.*(skperp2./ktot2).^2).*...
%         exp(-4*skperp.^2/sktheta.^2.*(sin(beta1/2)).^2);
%P_Bonoli1=((uzz+cos(beta1)*urr).^2+(sin(beta1)*uii).^2).*exp(-(lenls-sktheta).^2/sktheta^2);
% P_ott1(i,iii)=(1+2*skperp2./ktot2.*((wpe./wce).^2-(wpi./w).^2).*(sin(beta1(iii)/2)).^2+(sin(beta1(iii))).^2*wpe.^4/w/w/wce/wce.*(skperp2./ktot2).^2).*...
%         exp(-4*skperp.^2/sktheta.^2.*(sin(beta1(iii)/2)).^2);
    P_ott1_term1(i,iii)=(1+2*skperp2./ktot2.*((wpe./wce).^2-(wpi./w).^2).*(sin(beta1(iii)/2)).^2).*expbeta(i,iii);
        
    P_ott1_term2(i,iii)=(sin(beta1(iii))).^2*wpe.^4/w/w/wce/wce.*(skperp2./ktot2).^2.*expbeta(i,iii);
    P_ott1(i,iii)=P_ott1_term1(i,iii)+P_ott1_term2(i,iii);

%exp12=exp(-(skperp^2+skperm^2)/sktheta^2); 

%exp3=exp(-lenls2/sktheta^2);
%Vkkm2=0.25*w^2*abs((Hh*Hhm*cj)^-1);     %  |V|的求法
Vkkm2=0.25*w^2*abs((Hh*Hhm)^-1);     %  |V|的求法
Vkkm3(i,iii)=Vkkm2;
const=2*pi.*skperm./Vgperm.*Vkkm2.*(pi.*sktheta^2)^-1.*deltan^2;
const1(i,iii)=const;
%pks=const.*cj1.*exp3;
 P_Bonoli(i,iii)=((uzz+cos(beta1(iii))*urr).^2+(sin(beta1(iii))*uii).^2).*expbeta(i,iii);
 pkss(i,iii)=P_Bonoli(i,iii)*const;
 %P_Bonoli(i,iii)=const1(i,iii).*P_Bonoli1(i,iii);
 const_ott(i,iii)=2*pi.*skperm./Vgperm./deltaw(i).^2.*(pi.*sktheta^2)^-1.*deltan^2;
P_ott11(i,iii)=const_ott(i,iii).*P_ott1(i,iii);
P_ott11_term1(i,iii)=const_ott(i,iii).*P_ott1_term1(i,iii);
% P_ott=(1-eps.perp.*(sin(beta1/2)).^2+(sin(beta1)).^2*wpe.^4/w/w/wce/wce).*...
%         exp(-4*skperp.^2/sktheta.^2.*(sin(beta1/2)).^2);

%P_Bonoli(i,iii)=((uzz+cos(beta1(iii))*urr).^2+(sin(beta1(iii))*uii).^2).*exp(-(skperp^2+skperm^2-2*skperp*skperm*cos(beta1(iii)))/sktheta^2);
%b(i,iii)=P_Bonoli(i,iii)./const1(i,iii);
% nu=2*((sin(beta1(iii)/2)).^2*const.*P_Bonoli(i,iii));
nu=2*((sin(beta1(iii)/2)).^2.*P_ott11(i,iii));
% 
 nu1(i,iii)=nu;

%pkst1=@(beta2)(const./cj2*exp(-(skperp^2+skperm^2-2*skperp*skperm*cos(beta2))/sktheta^2)); %P（beta）
% pkst1=@(beta2)(const.*((uzz+cos(beta2)*urr).^2+(sin(beta2)*uii).^2).*exp(-(skperp^2+skperm^2-2*skperp*skperm*cos(beta2))/sktheta^2)); %P（beta）
% pkst=const.*cj1.*exp(-(skperp^2+skperm^2-2*skperp*skperm*cos(beta1))/sktheta^2);
% pkstot=integral(pkst1,-pi,pi);        %P（beta）对beta在（-pi，pi）上的积分
end
s_ott(i)=trapz(beta1,P_ott1(i,:));


    s(i)=trapz(beta1,P_Bonoli(i,:).*const1(i,:));
    s1(i)=trapz(beta1,nu1(i,:));
    ls(i,:)=Vgperm1(i,:)./s1(i);
     ratk=(skperp/sktheta)^2;
     if(ratk<1)
        pfs=ls(i,:)*sktheta;
    else
        pfs=ls(i,:)*sktheta/ratk;
     end
    pfs1(i,:)=pfs';

end
color_a = [0 0 1; 1 0 0;  1 0 1];

for i_len=1:length(npar1)
    if npar1(i_len)==-2.23 
            figure (111)
    pos=[0.15 0.1 0.25 0.4];
    subplot('Position',pos) 
        hold on;
        box on;
        h=plot((beta1)*180/pi,P_ott1(i_len,:),'-','linewidth',2,'color',color_a(color_ind,:));
        
        ylabel('P_{scatter} (normalized)');
        set(gca,'FontSize',20);
        set(findall(gcf,'type','line'),'linewidth',2);
        axis([-220 220 0 1.2])
        set(gca,'Xtick',[-180:90:180])
        set(gca,'Ytick',[0:0.6;1.2])
        set(gca,'ticklength',[0.015 0.015]);
         LA=axis;
        textstr={['total']};
        text(LA(1)+0.05*(LA(2)-LA(1)),LA(4)+0.05*(LA(3)-LA(4)),textstr,'FontSize',20);
        textstr2={['(d)']}
        text(LA(2)-0.1*(LA(2)-LA(1)),LA(4)+0.05*(LA(3)-LA(4)),textstr2,'FontSize',20);

        pos=[0.4 0.1 0.25 0.4];
    subplot('Position',pos) 
        hold on;
        box on;
        h=plot((beta1)*180/pi,P_ott1_term1(i_len,:),'-','linewidth',2,'color',color_a(color_ind,:));
        set(gca,'FontSize',20);
        set(findall(gcf,'type','line'),'linewidth',2);
        xlabel('{\beta} ( ^{o})');
        axis([-220 220 0 1.2])
        set(gca,'Xtick',[-180:90:180])
        set(gca,'Ytick',[0,0.6,1.2], 'YtickLabel',[]);
        set(gca,'ticklength',[0.015 0.015]);
                 LA=axis;
        textstr={['geometrical optics approximation']};
        text(LA(1)+0.05*(LA(2)-LA(1)),LA(4)+0.05*(LA(3)-LA(4)),textstr,'FontSize',20);

        textstr2={['(e)']}
        text(LA(2)-0.1*(LA(2)-LA(1)),LA(4)+0.05*(LA(3)-LA(4)),textstr2,'FontSize',20);

         pos=[0.65 0.1 0.25 0.4];
    subplot('Position',pos) 
        hold on;
        box on;       
        h=plot((beta1)*180/pi,P_ott1(i_len,:)-P_ott1_term1(i_len,:),'-','linewidth',2,'color',color_a(color_ind,:));
        set(gca,'FontSize',20);
        set(findall(gcf,'type','line'),'linewidth',2);
        axis([-220 220 0 1.2])
        set(gca,'Xtick',[-180:90:180])
        set(gca,'Ytick',[0, 0.6, 1.2], 'YtickLabel',[]);
        set(gca,'ticklength',[0.015 0.015]);
                 LA=axis;
        textstr={['E x B drifts']};
        text(LA(1)+0.05*(LA(2)-LA(1)),LA(4)+0.05*(LA(3)-LA(4)),textstr,'FontSize',20);
                text(LA(1)+0.05*(LA(2)-LA(1)),LA(4)+0.05*(LA(3)-LA(4)),textstr,'FontSize',20);
        textstr2={['(f)']}
        text(LA(2)-0.1*(LA(2)-LA(1)),LA(4)+0.05*(LA(3)-LA(4)),textstr2,'FontSize',20);
        
         pos=[0.15 0.5 0.25 0.4];
    subplot('Position',pos) 
        hold on;
        box on;
        h=plot((beta1)*180/pi,P_ott11(i_len,:)*tstep,'-','linewidth',2,'color',color_a(color_ind,:));
        
        ylabel('P_{scatter}');
        set(gca,'FontSize',20);
        set(findall(gcf,'type','line'),'linewidth',2);
        axis([-220 220 0 0.015])
        set(gca,'Xtick',[-180:90:180],'XtickLabel',[])
        set(gca,'Ytick',[0.005:0.005:0.015]);
        set(gca,'ticklength',[0.015 0.015]);
        set(0,'defaultfigurecolor','w'); 
        LA=axis;
        textstr={['total']};
        text(LA(1)+0.05*(LA(2)-LA(1)),LA(4)+0.05*(LA(3)-LA(4)),textstr,'FontSize',20);
        textstr2={['(a)']}
        text(LA(2)-0.1*(LA(2)-LA(1)),LA(4)+0.05*(LA(3)-LA(4)),textstr2,'FontSize',20);
        
     pos=[0.4 0.5 0.25 0.4];
    subplot('Position',pos) 
        hold on;
        box on;
        
        h=plot((beta1)*180/pi,P_ott11_term1(i_len,:)*tstep,'-','linewidth',2,'color',color_a(color_ind,:));

        set(gca,'FontSize',20);
        set(findall(gcf,'type','line'),'linewidth',2);
        axis([-220 220 0 0.015])
        set(gca,'Xtick',[-180:90:180],'XtickLabel',[])
        set(gca,'Ytick',[0.005:0.005:0.015], 'YtickLabel',[]);
        set(gca,'ticklength',[0.015 0.015]);
        set(0,'defaultfigurecolor','w'); 
        LA=axis;
        textstr={['geometrical optics approximation']};
        text(LA(1)+0.05*(LA(2)-LA(1)),LA(4)+0.05*(LA(3)-LA(4)),textstr,'FontSize',20);

        textstr2={['(b)']}
        text(LA(2)-0.1*(LA(2)-LA(1)),LA(4)+0.05*(LA(3)-LA(4)),textstr2,'FontSize',20);
        
          pos=[0.65 0.5 0.25 0.4];
    subplot('Position',pos) 
        hold on;
        box on;
        
        
        h=plot((beta1)*180/pi,(P_ott11(i_len,:)-P_ott11_term1(i_len,:))*tstep,'-','linewidth',2,'color',color_a(color_ind,:));

        set(gca,'FontSize',20);
        set(findall(gcf,'type','line'),'linewidth',2);
        axis([-220 220 0 0.015])
        set(gca,'Xtick',[-180:90:180],'XtickLabel',[])
        set(gca,'Ytick',[0.005:0.005:0.015], 'YtickLabel',[]);
        set(gca,'ticklength',[0.015 0.015]);
        set(0,'defaultfigurecolor','w'); 
        LA=axis;
        textstr={['E x B drifts']};
        text(LA(1)+0.05*(LA(2)-LA(1)),LA(4)+0.05*(LA(3)-LA(4)),textstr,'FontSize',20);
                text(LA(1)+0.05*(LA(2)-LA(1)),LA(4)+0.05*(LA(3)-LA(4)),textstr,'FontSize',20);
        textstr2={['(c)']};
        text(LA(2)-0.1*(LA(2)-LA(1)),LA(4)+0.05*(LA(3)-LA(4)),textstr2,'FontSize',20');

    end
end

%plot((beta1)*180/pi,P_ott11(i,:)./max(P_ott11(i,:)),'linewidth',2);
% legend_str{i}=['N_{//}=',num2str(npar)];
% legend(legend_str)
figure (999)
hold on
box on
plot(npar1,s./1e11,'-','linewidth',2,'color',color_a(color_ind,:))
xlabel('N_{//}');
set(gca,'Xlim',[-3 -1.6],'Xtick',[-2.84 -2.23 -1.82])
ylabel('Probability for {\beta}');
set(gca,'FontSize',20);
set(findall(gcf,'type','line'),'linewidth',2);
set(0,'defaultfigurecolor','w'); 


figure (9999)
hold on
plot(npar1,pfs1(:,1))



disp(s1./1e11)