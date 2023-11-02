ccc

r=4;
D=2;
c=sqrt(4*D*r)+1;
K=2;
x=linspace(0,K);
tspan = [0 100];
y0 = [K-0.1 -0.1];
[t,y] = ode45(@(t,y) odefcn(t,y,c,r,D,K), tspan, y0);

hold on
plot(x,-r/c*x.*(1-x/K),'k')
plot(x,-r/D*x,'k')
plot(y(:,1),y(:,2),'b')
axis([0 K min(-r/c*x.*(1-x/K)) 0])



function dydt = odefcn(t,y,c,r,D,K)

dydt = [y(2);
        -c*y(2)/D-r*y(1)*(1-y(1)/K)/D];
end