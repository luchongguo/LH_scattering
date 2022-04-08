fprintf('�����������½�:\n')
a=input('');
fprintf('�����������Ͻ�:\n')
b= input('');
fprintf('�������ֵalpha:\n')
alpha=input('');
fprintf('����������������N:\n')
N = input('');

x0=a;
y0=alpha;
h=(b-a)/N;
K=[];
x1=[];
y1=[];
for n=1:N
    K1=h*fun_Kutta(x0,y0);
    K2=h*fun_Kutta(x0+h/2,y0+K1/2);
    K3=h*fun_Kutta(x0+h/2,y0+K2/2);
    K4=h*fun_Kutta(x0+h,y0+K3/2);
    x1(n)=x0+h;
    y1(n)=y0+(K1+2*K2+2*K3+K4)/6;
    x0=x1(n);
    y0=y1(n);
end

function z=fun_Kutta(x,y)
 z = x + y;
% z = - y^2;
% z=2*y/x+x^2*exp(x);
% z=(y^2+y)/x;
% z=-20*(y-x^2)+2*x;
% z=-20*y+20*sin(x)+cos(x);
%z=-20*(y-exp(x)*sin(x))+exp(x)*(sin(x)+cos(x));
end
