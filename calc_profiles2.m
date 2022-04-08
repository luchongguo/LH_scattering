function [nete1]=calc_profiles2(ne,Ne,rho,rhopsi,rho1)

% [~,bb]=find(Ne==max(Ne));
% if Ne(bb)>ne(end)
%     ne(end)=Ne(bb);
%     ne(end-5:end)=interp1(rho(end-5:5:end),ne(end-5:5:end),rho(end-5:end));
% end
%     
% % find(max(Ne) )
% Ne=ne1;
% rhopsi=rhopsi11;
%  rho111=[0:0.001:1];
%  ne11=interp1(rho,ne,rho111);
%  [B,A]=min(abs(rhopsi(641,641:end)-1))
%
  rho111=[0:0.001:1];
  ne11=interp1(rho,ne,rho111);
%  [B,A]=min(abs(rhopsi(641,641:end)-1))
% 
% if B<rhopsi(641,640+A)    
%     ne11(A)=0;
% else
%     ne11(A)=0;
% end



  
 for i=1:1281
     
     [k,k1]=min(abs(rhopsi(641,i)-rho111));
        if rhopsi(641,i)<1
     nee(i)=ne11(k1);
        else 
            nee(i)=0;
    
        end
 end
nee2=nee+Ne;



for i=1:1281
    for j=1:1281
%         for k=1:101
            
                
                [~,k1]=min(abs(rhopsi(i,j)-rhopsi(641,:))); 
                nete1(i,j)=nee2(k1); 

            
%             if inpolygon(eqdsk_r(i),eqdsk_z(j),gvar.rlim,gvar.zlim) && inpolygon(eqdsk_r(i),eqdsk_z(j),gvar.rbbbs,gvar.zbbbs)~=0
%                 nete1(i,j)=0;
%             end 
    end
end