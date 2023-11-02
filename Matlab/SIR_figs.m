ccc
fs=15;
r=1;
a=2;

figure('position',[0 0 1/3 1/2])
hold on
plot(a/r*[1, 1],[0,5],'r--')
plot([0,5],[0,0],'b--')
plot([0,5],[5,0],'g--')

text(3,2,'$\dot{I}>0$','interpreter','latex','fontsize',fs)
text(3,3,'$\dot{S}<0$','interpreter','latex','fontsize',fs)
text(.5,2,'$\dot{I}<0$','interpreter','latex','fontsize',fs)
text(.5,3,'$\dot{S}<0$','interpreter','latex','fontsize',fs)
annotation('textarrow',[.55,.5]+0.1,[.55,.6]+0.2);
annotation('textarrow',[.55,.5]-0.15,[.6,.55]+0.2);

axis equal
axis([0 5 0 5])
xlabel('$S$ in thousands','interpreter','latex')
ylabel('$I$ in thousands','interpreter','latex')
legend('$I$ nullcline','$S$ and $I$ nullcline','$S_0+I_0=N$','location','s')
set(gca,'fontsize',fs)
export_fig('../Pictures/SIR_phase_plane.png','-r300')
%%
figure('position',[0 0 1/2 1/2])
hold on
plot(a/r*[1, 1],[0,5],'r--')
plot([0,5],[0,0],'b--')
plot([0,5],[5,0],'g--')
plot(a/r*[1, 1],[0,5],'--')
plot([0,5],[0,0],'--')
for i=1:4
[t1,y1]=ode45(@(t,y)ODE(t,y,r,a),[0,-100],[i 5-i 0]);
[t2,y2]=ode45(@(t,y)ODE(t,y,r,a),[0,100],[i 5-i 0]);
plot(y1(:,1),y1(:,2),':k')
plot(y2(:,1),y2(:,2),'-k')
end
axis equal
axis([0 10 0 5])
xlabel('$S$ in thousands','interpreter','latex')
ylabel('$I$ in thousands','interpreter','latex')
set(gca,'fontsize',fs)
export_fig('../Pictures/SIR_phase_plane_trajectories.png','-r300')
%%
figure('position',[0 0 1 1/2])
for i=1:4
    subplot(1,4,i)
[t2,y2]=ode45(@(t,y)ODE(t,y,r,a),[0,5],[i 5-i 0]);
plot(t2,y2)
axis([0 5 0 5])
xlabel('Time in days','interpreter','latex')
ylabel('Population in thousands','interpreter','latex')
legend('$S$','$I$','$R$','location','e')
set(gca,'fontsize',fs)
end
export_fig('../Pictures/SIR_trajectories.png','-r300')

function dydt=ODE(t,y,r,a)
dydt=[-r*y(1)*y(2);
       r*y(1)*y(2)-a*y(2);
       a*y(2)];
end