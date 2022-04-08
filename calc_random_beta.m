%calculate the random beta probability 
function [beta1]=calc_random_beta(uzz,urr,uii,skperp,skperm,ktheta)
%f1=@(beta2) 1./(1.2*pi*(1+5*(beta2-7.3).^2));
f1=@(beta2)(((uzz+cos(beta2)*urr).^2+(sin(beta2)*uii).^2).*...
    exp(-(skperp^2+skperm^2-2*skperp*skperm*cos(beta2))/ktheta^2)); %P��beta��
beta22=-pi:0.01:pi;
ff=f1(beta22);%���ݹ�ʽ��������ܶ�
s=trapz(beta22,ff);  %����������������ܶȵĻ���
%ff=ff/s;         %��һ�������ܶ�
 
 
while 1
    beta2=rand(1)*2*pi-pi;%����[0,24]���ȷֲ������
    if beta2<=pi
        f=f1(beta2)/max(ff);

    end         %�����Ӧ�ܶȺ���ֵf(beta2)
    r=rand(1);  %����[0,1]���ȷֲ������
    if r<=f     %��������rС��f(beta2)�����ɸ�t����������a��
        beta1=beta2;
        return
    end
end

%  cosbeta=(skperp^2+skperm^2-ktheta^2)./(2*skperp*skperm);
%  beta222=acos(cosbeta)
%  beta1=beta1+beta222;