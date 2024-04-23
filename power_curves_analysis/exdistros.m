% Check distros

figure, normplot(randn(1000,1))
figure, normplot(rand(1000,1))
figure, normplot(exprnd(1,1000,1))
x=randn(1000,1);
x(37)=10;
figure, normplot(x)