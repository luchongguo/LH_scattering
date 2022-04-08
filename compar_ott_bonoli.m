P_ott1=(1+2*skperp2./ktot2.*((wpe./wce).^2-(wpi./w).^2).*(sin(beta1/2)).^2+(sin(beta1)).^2*wpe.^4/w/w/wce/wce).*...
        exp(-4*skperp.^2/ktheta.^2.*(sin(beta1/2)).^2);
term1_ott=((1+2*skperp2./ktot2.*((wpe./wce).^2-(wpi./w).^2).*(sin(beta1/2)).^2)).^2;

term1_Bonoli=(uzz+cos(beta1)*urr).^2;

hold on;
plot(beta1,term1_ott./max(term1_ott));

plot(beta1,term1_Bonoli./max(term1_Bonoli));


a=term1_Bonoli./term1_ott;
figure
term2_ott=(sin(beta1)).^2*wpe.^4/w/w/wce/wce.*(skperp2./ktot2).^2;
term2_Bonoli=((sin(beta1)*uii).^2);
av=term2_Bonoli./term2_ott;
hold on
plot(beta1,term2_ott./max(term2_ott));
plot(beta1,term2_Bonoli./max(term2_Bonoli));