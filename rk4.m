function x=rk4(f,x0,h,a,b,N)
%x0��ֵ
%(a,b) Ϊ��������
%N��������
t=a:h:b;
%t=zeros(1:N+1);
x(:,1)=x0;
for i=1:N+1
L1=f(t(i),x(:,i));
L2=f(t(i)+h/2,x(:,i)'+(h/2)*L1);
L3=f(t(i)+h/2,x(:,i)'+(h/2)*L2);
L4=f(t(i)+h,x(:,i)'+h*L3);
x(:,i+1)=x(:,i)'+(h/6)*(L1+2*L2+2*L3+L4);
end

% f=@(t,x)[2*x(1)-0.02*x(1)*x(2),0.0002*x(1)*x(2)-0.8*x(2)];
% x0=[3000,120];
% N=10;
% a=0;b=10;
% rk4(f,x0,a,b,10)