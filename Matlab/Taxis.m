ccc
m = 0;
x = linspace(0,1,1000);
t = [0, 0.4, 0.8, 1];

sol = pdepe(m,@pdex1pde,@pdex1ic,@pdex1bc,x,t,odeset('RelTol',1e-7,'absTol',1e-7));
% Extract the first solution component as u.
u = sol(:,:,1);
% %%
% close all
% 
% for i=1:length(t)
%     hold off
%     plot(x,u(i,:))
%     hold on
%     plot(x,x)
%     
%     axis([0 1 0 2])
%     drawnow
% end
% 
% % A solution profile can also be illuminating.
% % figure
%%
close all
hold on
plot(x,u(1,:))
plot(x,u(2,:))
plot(x,u(3,:))
plot(x,u(4,:))
% plot(x,u(5,:))
plot(x,x)
xlabel('Distance, $x$')
ylabel('$u$')
format shortg;
legendCell = cellstr(num2str(round(t,2)', '$t$=%-.1f'))
legendCell{end+1,1}='$D_l=x$';
legend(legendCell)
set(gca,'fontsize',15)
xlim([0 1])
export_fig('../Pictures/Taxis.png','-r300')
% --------------------------------------------------------------
function [c,f,s] = pdex1pde(x,t,u,DuDx)
c = 1;
f = x*u;
s = 0;
end
% --------------------------------------------------------------
function u0 = pdex1ic(x)
u0 = exp(-((x-0.5)*10)^2);
end
% --------------------------------------------------------------
function [pl,ql,pr,qr] = pdex1bc(xl,ul,xr,ur,t)
pl = 0;
ql = 1;
pr = 0;
qr = 1;
end