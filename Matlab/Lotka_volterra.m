function Lotka_volterra(varargin)
% (c) Thomas E. Woolley, 2017.

clc
close all
fs=15;
alpha=1/2;
tspan = linspace(0,9,1e3); % Length of time the equation is solved for
figure('position',[0 1/3 2/3 1/3])
subplot(1,2,1)

for r=linspace(1,3,10)
u0=1;
v0=r;
[t,u] = ode23(@(t,u)ODE(t,u,alpha),tspan, [u0;v0]); % Matlab ODE solver.

hold on
plot(u(1,1),u(1,2),'ko','linewidth',3)
plot(u(:,1),u(:,2),'linewidth',3) % Plot solution.

end
axis([0 5 0 4])
xlabel('Hare','Interpreter','latex')
ylabel('Lynx','Interpreter','latex')
% L=legend('Hare','Lynx');
% set(L,'Interpreter','latex')
set(gca,'fontsize',fs)
subplot(1,2,2)

[t,u] = ode23(@(t,u)ODE(t,u,alpha),linspace(0,60,1e3), [1;3]); % Matlab ODE solver.

plot(t,u(:,1),'linewidth',3)
hold on
plot(t,u(:,2),'linewidth',3)
xlabel('Time','Interpreter','latex')
ylabel('Population','Interpreter','latex')
L=legend('Hare','Lynx');
set(L,'Interpreter','latex')
set(gca,'fontsize',fs)
set(0,'DefaultFigureColor',[1 1 1])
export_fig '../Pictures/LV_sims.png' -m3





function dydt = ODE(t,u,alpha)
% ODE equation to solve.
dydt = [u(1)-u(1)*u(2);
    alpha*u(2)*(u(1)-1)];
