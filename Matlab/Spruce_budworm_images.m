ccc
k=8;
x=linspace(0,k);
figure('position',[0 0 1 2/3])
l=0
for r=[0.25 0.5 0.7]
    l=l+1;
    subplot(2,3,l)
    hold on
    plot(x,r*x.*(1-x/k),'b')
    plot(x,x.^2./(1+x.^2),'r')
    axis([0 max(x) 0 1.5])
    
    L=legend('$ru(1-u/k)$','$u^2/(1+u^2)$');
    set(L,'interpreter','latex','location','south')
    xlabel('$u$','interpreter','latex')
    xticks([])
    yticks([])
    set(gca,'fontsize',15)
end

for r=[0.25 0.5 0.7]
    l=l+1;
    subplot(2,3,l)
    hold on
    %     plot(x,0*x,'k','linewidth',1)
    axis([0 max(x) -0.5 0.5])
    xlabel('$u$','interpreter','latex')
    ylabel('$\dot{u}$','interpreter','latex')
    PlotAxisAtOrigin(x,r*x.*(1-x/k)-x.^2./(1+x.^2))
    
    
    
    
    set(gca,'fontsize',15)
end