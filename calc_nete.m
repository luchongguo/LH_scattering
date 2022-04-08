function [nete1]=calc_NeTe(nete,rho,rhopsi)
for i=1:129
    for j=1:129
%         for k=1:101
            if rhopsi(i,j)>1
                nete1(i,j)=0;
            else
                 [~,k]=min(abs(rhopsi(i,j)-rho));
                 %[~,ii]=min(abs(z-eqdsk_z1));
                 nete1(i,j)=nete(k);
            end  
%         end  
    end
end