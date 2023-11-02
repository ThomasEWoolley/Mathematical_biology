ccc
fs=15;
options=odeset('reltol',1e-7);
for beta=[0.25 0.75 1 1.2 3]
e=1;
if beta==0.75
    e=0.01;
end
figure('position',[0 0 1/2 .25])
[t,y]=ode45(@(t,x)ode(t,x,beta),[0,30],[beta 1/beta+e],options);
subplot(1,2,1)
plot(t,y)
l=legend('$u$','$v$','location','best');
set(l,'interpreter','latex')
xlabel('$t$','interpreter','latex')
ylabel('$u$ or $v$','interpreter','latex')
set(gca,'fontsize',fs)
subplot(1,2,2)
hold on
plot(y(:,1),y(:,2))
plot(beta,1/beta,'o')
xlabel('$u$','interpreter','latex')
ylabel('$v$','interpreter','latex')
set(gca,'fontsize',fs)
% close all
export_fig(['..\Pictures\Schnakenberg_beta_',num2str(beta*100),'.png'],'-r300')
end


function dy=ode(t,x,beta)
u=x(1);
v=x(2);
dy=[-u+u^2*v;
    beta-u^2*v];
end