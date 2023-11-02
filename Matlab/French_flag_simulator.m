ccc
fs=20;
m = 0;n=10;
x = linspace(0,5,500);
t = linspace(0,20,n*20);


D=1;
gamma=0.1;
S=2;

u = pdepe(m,@(x,t,u,DuDx)pdex1pde(x,t,u,DuDx,D,gamma)...
    ,@pdex1ic,@(xl,ul,xr,ur,t)pdex1bc(xl,ul,xr,ur,t,S),x,t);
%%
close all
figure('position',[0 0.1 1 3/4])
subplot(1,2,1)
hold on
plot(x,u(1,:))
plot(x,u(n,:))
plot(x,u(2*n,:))
plot(x,u(5*n,:))
plot(x,u(10*n,:))
plot(x,u(20*n,:))
plot(x,S*cosh(sqrt(gamma/D)*(x-5))./cosh(sqrt(gamma/D)*(-5)),'g--')
xlabel('Distance $x$')
ylabel('$u$')
legend('$t=0$','$t=1$','$t=2$','$t=5$','$t=10$','$t=20$','Steady state solution','location','best')
set(gca,'fontsize',fs)
subplot(1,2,2)
hold on
surf(x,t,u)
xlabel('Distance, $x$')
ylabel('Time, $t$')
zlabel('$u$')
shading interp
axis tight
set(gca,'fontsize',fs)
export_fig('..\Pictures\French_flag_sims.png','-r300')
% Extract the first solution component as u.
%%
% filename='../Animations/Fisher.gif';
% if exist(filename, 'file')==2
%     delete(filename);
% end
% close all
% figure('position',[0 0.1 1 1/3])
% for i=1:length(t)
% plot(x,u(i,:))
% xlabel('Distance $x$')
% ylabel('$u$')
% set(gca,'fontsize',fs)
% drawnow
% frame = getframe(1);    im = frame2im(frame);    [imind,cm] = rgb2ind(im,256);
%     if i == 1;
%         imwrite(imind,cm,filename,'gif','Delaytime',0.1, 'Loopcount',inf);
%         for j=1:10
%             imwrite(imind,cm,filename,'gif','Delaytime',0.1,'WriteMode','append');
%         end
%     else
%         imwrite(imind,cm,filename,'gif','Delaytime',0.1,'WriteMode','append');
%     end
%      
% end
% for i=1:10
%     imwrite(imind,cm,filename,'gif','Delaytime',0.1,'WriteMode','append');
% end
% --------------------------------------------------------------
function [c,f,s] = pdex1pde(x,t,u,DuDx,D,gamma)
c = 1;
f = D*DuDx;
s = -gamma*u;
end
% --------------------------------------------------------------
function u0 = pdex1ic(x)
u0 = 0;
end
% --------------------------------------------------------------
function [pl,ql,pr,qr] = pdex1bc(xl,ul,xr,ur,t,S)
pl = ul-S;
ql = 0;
pr = 0;
qr = 1;
end