function [Ne,Te,delta,rho1,ktheta]=calc_edge(nedr,eqdsk_r,rhopsi,ne,te)
%����rhopsi�ķֲ�����Ne�ķֲ�
Ne=[];
[~,ii_start]=min(abs(eqdsk_r(641:end)-nedr(1,1)));
i=ii_start+640;
[~,iii_start]=min(abs(rhopsi(641,641:end)-1));
iii=iii_start+640;
 for ii=1:length(eqdsk_r)
%      [~,j]=min(abs(rhopsi(65,:)-1))
   % for jj=1:55
    %    if(abs(eqdsk_r(ii)-nedr(jj,1))<0.001)
    if eqdsk_r(ii)>nedr(1,1) & eqdsk_r(ii)<nedr(end,1)
        [~,k]=min(abs(eqdsk_r(ii)-nedr(:,1)));    
            Ne(ii)=nedr(k,2); 
            Te(ii)=nedr(k,3);
            delta(ii)=nedr(k,4);
            ktheta(ii)=nedr(k,5);
            rho1(ii)=rhopsi(641,ii);
             Ne(ii)=nedr(1,2)*exp(-(rho1(ii)-1)./0.05);
             %delta(ii)=nedr(1,4)*exp((rho1(ii)-1)./0.5);
             delta(ii)=nedr(1,4)*exp((rho1(ii)-1)./0.2);
%             
            
%             break;
    elseif ii>iii & ii<i 
        Ne(iii)=ne(end);
        
    else
            Ne(ii)=0;
            Te(ii)=0;
            rho1(ii)=0;
            delta(ii)=0;
            ktheta(ii)=0;
    end
%         
%     end
 end

 
 
