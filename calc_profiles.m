function [nete1]=calc_profiles(ne,Ne,rho,rhopsi,rho1)

% [~,bb]=find(Ne==max(Ne));
% if Ne(bb)>ne(end)
%     ne(end)=Ne(bb);
%     ne(end-5:end)=interp1(rho(end-5:5:end),ne(end-5:5:end),rho(end-5:end));
% end
%     
% % find(max(Ne) )
% Ne=ne1;
% rhopsi=rhopsi11;
 rho111=[0:0.001:1];
 ne1=interp1(rho,ne,rho111);

%  plot(rho111,ne1)
%  hold on
%  plot(rhopsi11(641,:),Ne)
%  [~,i]=max(Ne)
%  [~,j]=min(ne1)
%  rho12=[rho111(end-1) rhopsi11(641,i)]
%  [~,j]=min(abs(rhopsi(641,641:end)-1))
%  ne111(j+640:i)=interp1(rho12,[ne1(end-1) Ne(i)],rhopsi11(641,j+640:i))
%  plot(rho12,[ne1(end) Ne(i)])
%  plot()

for i=1:1281
    for j=1:1281
%         for k=1:101
            if rhopsi(i,j)>1&rhopsi(i,j)<1.8
                
                [~,k1]=min(abs(rhopsi(i,j)-rho1)); 
                nete1(i,j)=Ne(k1); 
                
            elseif rhopsi(i,j)<1 | rhopsi(i,j)==1
                 [~,k]=min(abs(rhopsi(i,j)-rho111));
                 %[~,ii]=min(abs(z-eqdsk_z1));
                 nete1(i,j)=ne1(k);
            else
                nete1(i,j)=0;
            end
%             if inpolygon(eqdsk_r(i),eqdsk_z(j),gvar.rlim,gvar.zlim) && inpolygon(eqdsk_r(i),eqdsk_z(j),gvar.rbbbs,gvar.zbbbs)~=0
%                 nete1(i,j)=0;
%             end 
    end
end