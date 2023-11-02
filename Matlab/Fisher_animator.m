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
filename='../Animations/Fisher.gif';
if exist(filename, 'file')==2
    delete(filename);
end
close all
figure('position',[0 0.1 1 1/3])
for i=1:length(t)
plot(x,u(i,:))
xlabel('Distance $x$')
ylabel('$u$')
set(gca,'fontsize',fs)
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