Nperp=0.5;
Npar=0.7;
X.Xxx=0;
X.Xxy=0;
X.Xxz=0;
X.Xyy=0;
X.Xyz=0;
X.Xzz=0;


Nperp1 = Nperp;
Nperp2 = Nperp1*Nperp1;    
Nperp3 = Nperp1*Nperp1*Nperp1;    
Nperp4 = Nperp1*Nperp1*Nperp1*Nperp1;
Npar1 = Npar;
Npar2 = Npar1*Npar1;
Npar3 = Npar1*Npar1*Npar1;
Npar4 = Npar1*Npar1*Npar1*Npar1;
N2=Nperp2+Npar2;

N6 = N2*N2*N2;

Xxx = X.Xxx;
Xxy = X.Xxy;
Xxz = X.Xxz;
Xyy = X.Xyy;
Xyz = X.Xyz;
Xzz = X.Xzz;
Xxx2 = X.Xxx*X.Xxx;
Xxy2 = X.Xxy*X.Xxy;
Xxz2 = X.Xxz*X.Xxz;
Xyy2 = X.Xyy*X.Xyy;
Xyz2 = X.Xyz*X.Xyz;
Xzz2 = X.Xzz*X.Xzz;

D(43) = Nperp4/N6;
D(1) = Xxx*Nperp4/N6;
D(2) = 2.0*Npar1*Xxz*Nperp3/N6;
D(3) = Npar2*Xxx*Nperp2/N6;
D(4) = Npar2*Xzz*Nperp2/N6;
D(5) = 2.0*Npar2*Nperp2/N6;
D(6) = -2.0*Xxx*Nperp2/N6;
D(7) = -Xyy*Nperp2/N6;
D(8) = -Xzz*Nperp2/N6;
D(9) = -Xxx*Xzz*Nperp2/N6;
D(10) = Xxz2*Nperp2/N6;
D(11) = -Xxx*Xyy*Nperp2/N6;
D(12) = -Xxy2*Nperp2/N6;
D(13) = -2.0*Nperp2/N6;    
D(14) = -2.0*Npar1*Xyy*Xxz*Nperp1;
D(15) = 2.0*Npar1*Xxy*Xyz*Nperp1;
D(16) = -2.0*Npar1*Xxz*Nperp1/N6;
D(17) = 2.0*Npar3*Xxz*Nperp1/N6;
D(18) = 1/N6;
D(19) = Npar4*Xzz/N6;
D(20) = Npar4/N6;
D(21) = Npar2*Xxz2/N6;
D(22) = -Npar2*Xyy*Xzz/N6;
D(23) = -Npar2*Xxx*Xzz/N6;
D(24) = -2.0*Npar2*Xzz/N6;
D(25) = -Npar2*Xyy/N6;
D(26) = -Npar2*Xyz2/N6;
D(27) = -Npar2*Xxx/N6;
D(28) = -2.0*Npar2/N6;
D(29) = Xyy*Xzz/N6;
D(30) = Xxx*Xzz/N6;
D(31) = Xxx*Xyy/N6;
D(32) = Xxx*Xyz2/N6;
D(33) = Xyz2/N6;
D(34) = Xxx/N6;
D(35) = Xxy2/N6;
D(36) = Xxy2*Xzz/N6;
D(37) = -Xxz2/N6;
D(38) = -Xyy*Xxz2/N6;
D(39) = Xyy/N6;
D(40) = Xzz/N6;
D(41) = Xxx*Xyy*Xzz/N6;
D(42) = 2.0*Xxy*Xyz*Xxz/N6;



