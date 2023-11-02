
clear all
close all
clc


fs=15; % Set fontsize

N=linspace(0,4);

for r=[1 2.1 3]
    figure
    hold on
    plot(N,N,'k') % Plot N(t)=N(t+1)
    plot(N,N.*exp(r*(1-N)),'b') % Plot N(t+1)=f(N(t))
    axis([0 4 0 3])
    xlabel('$N_t$')
    ylabel('$N_{t+1}$')
    
    Cobweb(0.1,r,'r')
    Cobweb(4,r,'g')
    Cobweb(3,r,'m')
    set(gca,'fontsize',fs) % Set fontsize.
    export_fig(['../Pictures/Ricker_cobweb_r',num2str(r),'.png'],'-r300')
end
function Cobweb(N0,r,c)
% The following for loop computes the first 100 iteration values.
Nt(1)=N0;
for i=1:100
    Nt(i+1)=Nt(i)*exp(r*(1-Nt(i)));
end

% The following code plots the cobweb.
plot([Nt(1) Nt(1)],[0 Nt(2)],c,'linewidth',1)
for i=1:100-2
    plot([Nt(i) Nt(i+1)],[Nt(i+1) Nt(i+1)],c,'linewidth',1)
    plot([Nt(i+1) Nt(i+1)],[Nt(i+1) Nt(i+2)],c,'linewidth',1)
end
end