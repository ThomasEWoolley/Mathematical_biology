ccc
m = 0;
x = linspace(0,1,1000);
t = [0 0.01 0.05 0.1 1];

sol = pdepe(m,@pdex1pde,@pdex1ic,@pdex1bc,x,t,odeset('RelTol',1e-7,'absTol',1e-7));
% Extract the first solution component as u.
u = sol(:,:,1);

% % A surface plot is often a good way to study a solution.
% surf(x,t,u()) 
% xlabel('Space, x')
% ylabel('Time, t')

% A solution profile can also be illuminating.
% figure
%%
close all
hold on
plot(x,u(1,:))
plot(x,u(2,:))
plot(x,u(3,:))
plot(x,u(4,:))
plot(x,u(5,:))
xlabel('Distance, $x$')
ylabel('$u$')
format shortg;
legendCell = cellstr(num2str(round(t,2)', '$t$=%-.2f'))
legend(legendCell)
set(gca,'fontsize',15)
xlim([0 1])
export_fig('../Pictures/Bounded_diffusion_zero_flux.png','-r300')
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
function [pl,ql,pr,qr] = pdex1bc(xl,ul,xr,ur,t)
pl = 0;
ql = 1;
pr = 0;
qr = 1;
end