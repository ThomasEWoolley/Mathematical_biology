ccc
cols='gbm';
fs=20;
T=10;
rvec=[0.9 1.5 3];
figure('position',[0 0 1 .75])
x0=0.5;
x=linspace(0,1);
for i=1:3
    subplot(2,3,i)
    hold on
%     plot(x,x,'k')
%     plot(x,rvec(i)*x.*(1-x),'b')
    f=[num2str(rvec(i)),'*x.*(1-x)'];
    cobweb(f, [0 1], x0, 1, 10,cols(i));
    xlabel('$N_t$')
    ylabel('$N_{t+1}$')
    set(gca,'fontsize',fs)
end



for i=1:3
    subplot(2,3,[4:6])
    hold on
    r=rvec(i);
    n(1)=x0;
    for j=1:T
        n(j+1)=r*n(j)*(1-n(j));
    end
    plot(0:T,n,'-o','color',cols(i));
    
    
end
xlabel('$t$')
ylabel('$N_t$')
set(gca,'fontsize',fs)
axis tight
legend('$r=0.9$','$r=1.5$','$r=3$','location','best');%[0.7 0.25 0.2 0.2])

% export_fig('../Pictures/Logistic_difference_equation.png','-r300')