a=0;
c=pi;
z=0;
r=rand(1,10);
parfor i=1:10;
    a=i;
    z=z+i;
    b(i)=r(i);
    if i<=c
        d=2*a;
    end
end