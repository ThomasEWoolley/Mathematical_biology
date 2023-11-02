ccc

t=0;
x=linspace(-100,100,1e6);
c=1/sqrt(2);
col=['r','b','k']
i=1;
for t=[-2 0 10]
u=1-exp(c*(x-c*t));
xp=x(u>0);
up=u(u>0);
xn=x(u<0);
un=u(u<0);
hold on
axis([-10 10 0 1])
plot(xp,up,col(i))
plot(xn,un*0,col(i))

i=i+1;
end
xlabel('$x$')
ylabel('$u$')
legend('$t=-2$','$t=0$','$t=10$')
export_fig('../Pictures/Exponential_travelling_wave.png','-r300')