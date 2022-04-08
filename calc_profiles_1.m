function [nete1]=calc_profiles(ne,Ne,rho,rhopsi,rho1)
% rho111=[0:0.001:1];
% ne1=interp1(rho,ne,rho111);
% rhopsi=rhopsi11;
% Ne=delta1;
[B,A]=min(abs(rhopsi(641,641:end)-1))

if B<rhopsi(641,640+A)    
    ne(end)=0;
else
    ne(end)=0;
end
 rho111=[0:0.001:1];
 ne11=interp1(rho,ne,rho111);
 
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

% 
% for i=1:1281
%     for j=1:1281
%    if rhopsi(i,j)>1&rhopsi(i,j)<1.8
%                 [~,k1]=min(abs(rhopsi(i,j)-rho1));
%                 nete1(i,j)=Ne(k1);
%                 
%             elseif rhopsi(i,j)<1 | rhopsi(i,j)==1
%                  [~,k]=min(abs(rhopsi(i,j)-rho111));
%                  %[~,ii]=min(abs(z-eqdsk_z1));
%                  nete1(i,j)=ne1(k);
%             else
%                 nete1(i,j)=0;
%             end
% %             if inpolygon(eqdsk_r(i),eqdsk_z(j),gvar.rlim,gvar.zlim) && inpolygon(eqdsk_r(i),eqdsk_z(j),gvar.rbbbs,gvar.zbbbs)~=0
% %                 nete1(i,j)=0;
% %             end 
%     end
% end