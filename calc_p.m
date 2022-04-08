
N=100000; %��Ҫ������ĸ���
a=zeros(N,1); %��������������
n=0;
%f1=@(beta2) 1./(1.2*pi*(1+5*(beta2-7.3).^2));
f1=@(beta2)(const.*((uzz+cos(beta2)*urr).^2+(sin(beta2)*uii).^2).*...
    exp(-(skperp^2+skperm^2-2*skperp*skperm*cos(beta2))/ktheta^2)); %P��beta��
beta22=-pi:0.01:pi;
ff=f1(beta22);%���ݹ�ʽ��������ܶ�
s=trapz(beta22,ff);  %����������������ܶȵĻ���
ff=ff/s;         %��һ�������ܶ�
 
 
while n<N
    beta2=rand(1)*2*pi-pi;%����[0,24]���ȷֲ������
    if beta2<=pi
        f=f1(beta2)/s;

    end         %�����Ӧ�ܶȺ���ֵf(beta2)
    r=rand(1);  %����[0,1]���ȷֲ������
    if r<=f     %��������rС��f(beta2)�����ɸ�t����������a��
        n=n+1;
        a(n)=beta2;
    end
end
 
%����Ϊ�����������a�Ĺ��̣�����Ϊͳ�Ƽ�����������Ƿ���Ϸֲ� 
num=100;         %��100������ͳ��
[x,c]=hist(a,num);    %ͳ�Ʋ�ͬ������ֵĸ���
dc=2*pi/num;        %�����С
x=x/N/dc;         %����ͳ�ƽ����������ܶ�
 
bar(c,x,1); hold on;  %����ͳ�ƽ���������ܶ�ֱ��ͼ
plot(beta22,ff,'r'); hold off; %���ݹ�ʽ�������ܶ�����
