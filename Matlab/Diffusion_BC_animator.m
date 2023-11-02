ccc
m = 0;
x1 = linspace(0,1,1000);
x2 = linspace(-10,10,2000);
t = linspace(0,0.4,500);

sol1 = pdepe(m,@pdex1pde,@pdex1ic,@pdex1bc1,x2,t,odeset('RelTol',1e-7,'absTol',1e-7));
sol2 = pdepe(m,@pdex1pde,@pdex1ic,@pdex1bc1,x1,t,odeset('RelTol',1e-7,'absTol',1e-7));
sol3 = pdepe(m,@pdex1pde,@pdex1ic,@pdex1bc2,x1,t,odeset('RelTol',1e-7,'absTol',1e-7));
% Extract the first solution component as u.
u1 = sol1(:,:,1);
u2 = sol2(:,:,1);
u3 = sol3(:,:,1);


%%
close all
filename='../Animations/Diffusions.gif';
if exist(filename, 'file')==2
    delete(filename);
end
figure('position',[0 0 1 1/3])
for i=1:length(t)
subplot(1,3,1)

plot(x2,u1(i,:))
title('Unbounded domain')
formatter
subplot(1,3,2)

plot(x1,u2(i,:))
text(.5,.9,['Time=',num2str(round(t(i),2))],'fontsize',15,'HorizontalAlignment','center')
title('Zero-flux conditions')
formatter
subplot(1,3,3)

plot(x1,u3(i,:))
title('Zero-Dirichlet conditions')
formatter
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

function formatter
xlabel('Distance, $x$')
ylabel('$u$')
set(gca,'fontsize',15)
axis([0 1 0 1])
end

% --------------------------------------------------------------
function [c,f,s] = pdex1pde(x,t,u,DuDx)
c = 1;
f = DuDx;
s = 0;
end
% --------------------------------------------------------------
function u0 = pdex1ic(x)
u0 = exp(-((x-1/2)*10)^2);
end
% --------------------------------------------------------------
function [pl,ql,pr,qr] = pdex1bc1(xl,ul,xr,ur,t)
pl = 0;
ql = 1;
pr = 0;
qr = 1;
end
function [pl,ql,pr,qr] = pdex1bc2(xl,ul,xr,ur,t)
pl = ul;
ql = 0;
pr = ur;
qr = 0;
end