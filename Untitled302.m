
x=1:129;
y=1:129;

[x,y]=meshgrid(x,y);
figure(1);
mesh(x,y,Bt);
xlabel('x');
ylabel('y');

figure(2);
xi=1:0.1:129;
yi=1:0.1:129;
[xi,yi]=meshgrid(xi,yi);
zi=interp2(x,y,Bt,xi,yi,'cubic');
mesh(xi,yi,zi);
xlabel('x');
ylabel('y');
