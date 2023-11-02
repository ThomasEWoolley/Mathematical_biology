ccc
fs=15;
r=1;
a=2;


figure('position',[0 0 1/2 1/2])
hold on
t=linspace(0,10,1e5);
for i=1:4
[t1,y1]=ode45(@(t,y)ODE(t,y,2,a),t,[4.9 0.1 0]);
[t2,y2]=ode45(@(t,y)ODE(t,y,1,a),t,[4.9 0.1 0]);
plot(t,y1(:,2),'-r')
plot(t,y2(:,2),'-b')
end
trapz(t,y1(:,2))
trapz(t,y2(:,2))
y1(end,3)/a
y2(end,3)/a
(5-y1(end,1))/a
(5-y2(end,1))/a
axis equal
axis([0 5 0 3])
xlabel('Time, $t$','interpreter','latex')
ylabel('$I$ in thousands','interpreter','latex')
legend('$r=2$','$r=1$')
set(gca,'fontsize',fs)

function dydt=ODE(t,y,r,a)
dydt=[-r*y(1)*y(2);
       r*y(1)*y(2)-a*y(2);
       a*y(2)];
end