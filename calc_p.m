
N=100000; %需要随机数的个数
a=zeros(N,1); %存放随机数的数列
n=0;
%f1=@(beta2) 1./(1.2*pi*(1+5*(beta2-7.3).^2));
f1=@(beta2)(const.*((uzz+cos(beta2)*urr).^2+(sin(beta2)*uii).^2).*...
    exp(-(skperp^2+skperm^2-2*skperp*skperm*cos(beta2))/ktheta^2)); %P（beta）
beta22=-pi:0.01:pi;
ff=f1(beta22);%根据公式计算概率密度
s=trapz(beta22,ff);  %计算整个区间概率密度的积分
ff=ff/s;         %归一化概率密度
 
 
while n<N
    beta2=rand(1)*2*pi-pi;%生成[0,24]均匀分布随机数
    if beta2<=pi
        f=f1(beta2)/s;

    end         %计算对应密度函数值f(beta2)
    r=rand(1);  %生成[0,1]均匀分布随机数
    if r<=f     %如果随机数r小于f(beta2)，接纳该t并加入序列a中
        n=n+1;
        a(n)=beta2;
    end
end
 
%以上为生成随机数列a的过程，以下为统计检验随机数列是否符合分布 
num=100;         %分100个区间统计
[x,c]=hist(a,num);    %统计不同区间出现的个数
dc=2*pi/num;        %区间大小
x=x/N/dc;         %根据统计结果计算概率密度
 
bar(c,x,1); hold on;  %根据统计结果画概率密度直方图
plot(beta22,ff,'r'); hold off; %根据公式画概率密度曲线
