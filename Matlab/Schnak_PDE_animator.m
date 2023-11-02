ccc
fs=15;
alpha=0.1;
beta=0.9;
Dv=20;
L=15;
x = linspace(0,L,10*L);
t = linspace(0,100,500);

sol = pdepe(0,@(x,t,u,DuDx)pdex1pde(x,t,u,DuDx,alpha,beta,Dv)...
    ,@(x)pdex1ic(x,alpha,beta),@pdex1bc,x,t);
% Extract the first solution component as u.
u = sol(:,:,1);
v = sol(:,:,2);
m1=min(min([u(:),v(:)]));
m2=max(max([u(:),v(:)]));
figure('position',[0 0 .5 1/3])
subplot(1,2,1)
pcolor(x,t,u)
shading interp
colorbar
xlabel('Distance, $x$')
ylabel('Time, $t$')
set(gca,'fontsize',fs)
subplot(1,2,2)
for i=1:length(t)
hold off
plot(x,u(i,:))
hold on
plot(x,v(i,:))
axis([0 L m1 m2])
title(['Time = ',num2str(round(t(i),1))])
xlabel('Distance, $x$')
legend('$u$','$v$')
set(gca,'fontsize',fs)
drawnow
end

% --------------------------------------------------------------
function [c,f,s] = pdex1pde(x,t,u,DuDx,alpha,beta,Dv)
c = [1;1];
f = [1;Dv].*DuDx;
s = [alpha-u(1)+u(1)^2*u(2);beta-alpha-u(1)^2*u(2)];
end
% --------------------------------------------------------------
function u0 = pdex1ic(x,alpha,beta)
u0 = [beta+rand*0.1;(beta-alpha)/beta^2+rand*0.1];
end
% --------------------------------------------------------------
function [pl,ql,pr,qr] = pdex1bc(xl,ul,xr,ur,t)
pl = [0;0];
ql = [1;1];
pr = [0;0];
qr = [1;1];
end