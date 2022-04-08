date=datetime('now');
clear global;
w11=2.45;
w111=num2str(w11);
ktheta111=num2str(ktheta11);
ktheta111=strrep(ktheta111,'.','p');
kmean1=num2str(kmean);
kmean1=strrep(kmean1,'.','p');
w111=strrep(w111,'.','p');
nnpar=num2str(npar);
nnpar=strrep(nnpar,'.','p');
zz1=num2str(zz);
zz1=strrep(zz1,'.','p');
lam_ne1=num2str(lam_ne);
lam_ne1=strrep(lam_ne1,'.','p');
lam_delta1=num2str(lam_delta);
lam_delta1=strrep(lam_delta1,'.','p');
ne_LCFS=ceil(max(ne1)./1e18);
tstep1=tstep./1e-11;


% savePath1 = strcat('/gpfs/scratch/chbwu/density_fluction_to_shenma/densityflucdata/90331/z_',num2str(zz1),'/picture');
% 
% mkdir(savePath1);
savePath2 = strcat('E:\OneDrive - mail.ustc.edu.cn\¸öÈËĞÅÏ¢\Ñ§¼®Ïà¹Ø\¸öÈËÂÛÎÄ\ÃÜ¶ÈÕÇÂä\LH_scattering_code\result\',num2str(nray),'_l',num2str(lengthray),'_tstep',num2str(tstep1));
mkdir(savePath2);
savePath3 = strcat('\S',num2str(shotno),'_r',num2str(rr),'_f',w111,...
    '_Npar',nnpar,'_ne',num2str(ne_LCFS),'e18_lam_ne',lam_ne1,'_delta',lam_delta1,'_kmean',num2str(kmean1),'_ktheta',num2str(ktheta111)) ;  
% savePath = strcat('F:\cbwu\ï¿½ì¹«\density_fluction_0408\densityflucdata\m',...
%     num2str(date.Month),'d',num2str(date.Day),'n',...
%     num2str(nray),'l',num2str(lengthray),'d_',...
%     w111,'f_',nnpar,'.mat'); % Æ´ï¿½ï¿½Â·ï¿½ï¿½ï¿½ï¿½ï¿½Ä¼ï¿½ï¿½ï¿½
save([savePath2,savePath3,'.mat']); % ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½val1,val2ï¿½ï¿½1_result.matï¿½ï¿½)

%disp(['delta=',num2str(delta)]);
disp(['nary=',num2str(nray)]);
disp(['lengthray=',num2str(lengthray)]); 
disp(['f=',num2str(w11)]);
disp(['Npar=',num2str(npar)]);
disp(['ktheta=',num2str(ktheta11)])
