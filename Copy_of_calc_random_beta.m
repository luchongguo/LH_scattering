%calculate the random beta probability 
function [beta1]=Copy_of_calc_random_beta(uzz,urr,uii,skperp,skperm,ktheta,kmean)

beta22=-pi:0.01:pi;
kmean2=kmean*200*pi;

lenls=sqrt(skperp^2+skperm^2-2*skperp*skperm*cos(beta22));
for iii=1:length(beta22)
    if beta22(iii)<0
        lenls(iii)=-abs(lenls(iii));
    else
        lenls(iii)=abs(lenls(iii));
    end
end

PP=((uzz+cos(beta22)*urr).^2+(sin(beta22)*uii).^2).*exp(-(lenls-kmean2).^2/ktheta^2);


%ff=f1(beta22);%���ݹ�ʽ��������ܶ�
s=trapz(beta22,PP);  %����������������ܶȵĻ���
%ff=ff/s;         %��һ�������ܶ�
 
 i=0;
while 1
    i=i+1;
    beta2=rand(1)*2*pi-pi;
    [~,ind]=min(abs(beta2-beta22));
    if beta2<=pi
        f=PP(ind)/max(PP);

    end         %�����Ӧ�ܶȺ���ֵf(beta2)
    r=rand(1);  %����[0,1]���ȷֲ������
    if r<=f     %��������rС��f(beta2)�����ɸ�t����������a��
        beta1=beta2;
        return  
    elseif i==1000
        beta1=0;
        return
    end
end






%  cosbeta=(skperp^2+skperm^2-ktheta^2)./(2*skperp*skperm);
%  beta222=acos(cosbeta)
%  beta1=beta1+beta222;