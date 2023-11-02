ccc
fs=20;
m = 0;n=10;
x = linspace(0,30,300);
t = linspace(0,10,n*10);

r=4;
D=1;
K=2;
c=sqrt(4*D*r);

u = pdepe(m,@(x,t,u,DuDx)pdex1pde(x,t,u,DuDx,r,D,K),@pdex1ic,@(xl,ul,xr,ur,t)pdex1bc(xl,ul,xr,ur,t,r,D,K),x,t);
% Extract the first solution component as u.
%%
close all
figure('position',[0 0.1 1 1/3])
subplot(1,3,1)
hold on
plot(x,u(1,:))
plot(x,u(n,:))
plot(x,u(2*n,:))
plot(x,u(5*n,:))
plot(x,u(10*n,:))
xlabel('Distance $x$')
ylabel('$u$')
legend('$t=0$','$t=1$','$t=2$','$t=5$','$t=10$','location','best')
set(gca,'fontsize',fs)

subplot(1,3,2)
hold on
v=linspace(0,K);
Times=t'*ones(1,length(x));
Space=ones(length(t),1)*x;
z=Space-c*Times;
dz=diff(z,1,2);
du=diff(u,1,2);
dudz=du./dz;
p3=plot(u(5*n,2:end),dudz(5*n,:),'g');
p1=plot(v,-r/c*v.*(1-v/K),'b--');
plot(v,0*v,'b--')
p2=plot(v,0*v,'r-','linewidth',1);
plot(v,-c/D*v,'r-','linewidth',1)
plot([K K],[0 -c*K/D],'r-','linewidth',1)
shading interp
xlabel('$u$')
ylabel('d$u$/d$z$')
axis([0 2 -1 0])
legend([p3,p1,p2],'Trajectory','Nullclines','$R$','location','best')
set(gca,'fontsize',fs)

subplot(1,3,3)
hold on
surf(t,x,u')
ylabel('Distance, $x$')
xlabel('Time, $t$')
zlabel('$u$')
shading interp
p1=plot3(t,c*t,3+0*t,'g:');
legend(p1,'Line of gradient $c$','location','se')
axis tight
set(gca,'fontsize',fs)
axis([0 10 0 30])
export_fig('..\Pictures\Fisher_sims.png','-r300')
% --------------------------------------------------------------
function [c,f,s] = pdex1pde(x,t,u,DuDx,r,D,K)
c = 1;
f = D*DuDx;
s = r*u*(1-u/K);
end
% --------------------------------------------------------------
function u0 = pdex1ic(x)
u0 = 0;
end
% --------------------------------------------------------------
function [pl,ql,pr,qr] = pdex1bc(xl,ul,xr,ur,t,r,D,K)
pl = ul-K;
ql = 0;
pr = ur;
qr = 0;
end