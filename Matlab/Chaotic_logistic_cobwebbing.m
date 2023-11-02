ccc
fs=20;
T=30;
rvec=[3.6];
figure('position',[0 0 1 .7])
x0=0.5;
x=linspace(0,1);

subplot(1,3,1)
hold on
%     plot(x,x,'k')
%     plot(x,rvec(i)*x.*(1-x),'b')
f=[num2str(rvec(1)),'*x.*(1-x)'];
cobweb(f, [0 1], x0, 1, 10,'b');
% cobweb(f, [0 1], x0+.1, 1, 10,'g');
xlabel('$N_t$')
ylabel('$N_{t+1}$')
set(gca,'fontsize',fs)
%%
subplot(1,3,[2:3])
hold on
r=rvec(1);
n(1)=x0;
for j=1:T
    n(j+1)=r*n(j)*(1-n(j));
end
plot(0:T,n,'-o','color','b');
% n(1)=x0+.1;
% for j=1:T
%     n(j+1)=r*n(j)*(1-n(j));
% end
% plot(0:T,n,'-o','color','g');
xlabel('$t$')
ylabel('$N_t$')
set(gca,'fontsize',fs)

% export_fig('../Pictures/Logistic_difference_equation_chaos.png','-r300')