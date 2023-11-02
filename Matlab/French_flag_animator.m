ccc
fs=20;
m = 0;n=10;
x = linspace(0,5,500);
t = linspace(0,20,1e3);


D=1;
gamma=0.1;
S=2;

u = pdepe(m,@(x,t,u,DuDx)pdex1pde(x,t,u,DuDx,D,gamma)...
    ,@pdex1ic,@(xl,ul,xr,ur,t)pdex1bc(xl,ul,xr,ur,t,S),x,t);

%%
filename='../Animations/French_flag.gif';
if exist(filename, 'file')==2
    delete(filename);
end
close all
figure('position',[0 0.1 1 1/3])
for i=1:length(t)
    
    hold off
plot(x,u(i,:))
hold on
plot(x,S*cosh(sqrt(gamma/D)*(x-5))./cosh(sqrt(gamma/D)*(-5)),'g--')
xlabel('Distance $x$')
ylabel('$u$')
set(gca,'fontsize',fs)
axis([0 5 0 2])
title(['$t=$',num2str(round(t(i),1))])
drawnow

frame = getframe(1);    im = frame2im(frame);    [imind,cm] = rgb2ind(im,256);
    if i == 1;
        imwrite(imind,cm,filename,'gif','Delaytime',0.1, 'Loopcount',inf);
        for j=1:10
            imwrite(imind,cm,filename,'gif','Delaytime',0.1,'WriteMode','append');
        end
    else
        imwrite(imind,cm,filename,'gif','Delaytime',0.1,'WriteMode','append');
    end
     
end
for i=1:10
    imwrite(imind,cm,filename,'gif','Delaytime',0.1,'WriteMode','append');
end
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