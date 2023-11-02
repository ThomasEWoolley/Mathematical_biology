ccc
fs=20;
a=0:0.001:15;
N=2000; x=0.1;
%%
close all
figure('position',[0 0 1 1/2])
hold on
for n=0:N
    x=a.^3.*x./(1+x).^3;
    
    if n>.75*N
        plot(a,x,'b.','MarkerSize',3)
    end
    
end
hold off
ylabel('Stable trajectories (log scale)')
xlabel('$\mu$')
xticks(1:15)
set(gca,'fontsize',fs)
set(gca,'yscale','log')
ylim([0.01,600])
export_fig('../Pictures/Discrete_inverse_cubic.png','-r300')
 