ccc
fs=20;
a=0:0.01:4;
N=100; x=0.1;

close all
figure('position',[0 0 1 1/2])
subplot(1,2,1)
hold on
a_vec=[];
x_vec=[];
for n=0:N
    x=a.*x.*(1-x);
    
    if n>.5*N
        plot(a,x,'b.','MarkerSize',3)
        a_vec=[a_vec,a];
        x_vec=[x_vec,x];
    end
    
end

%%
m = [a_vec;x_vec]; 
writematrix(m','C:\Users\smatw1\Documents\R\Discrete_logistic\Logistic_points.csv') 

plot([0 1],[0 0],'k')
% plot([1 4],[0 0],'--k')
x1=linspace(1,3);
x2=linspace(3,1+sqrt(6));
plot(x1,1-1./x1,'k')
plot(x2,(1+x2)./(2*x2).*(1+sqrt(1-4./(x2+1))),'r')
plot(x2,(1+x2)./(2*x2).*(1-sqrt(1-4./(x2+1))),'r')
axis ([0 4 0 1])
hold off
ylabel('$N^*$')
xlabel('$r$')
set(gca,'fontsize',fs)

subplot(1,2,2)
hold on
for n=0:N
    x=a.*x.*(1-x);
    if n>.5*N
        plot(a,x,'b.','MarkerSize',1)
    end
end
axis ([3.5 4 0 1])
hold off
ylabel('$N^*$')
xlabel('$r$')
set(gca,'fontsize',fs)
% export_fig('../Pictures/Logistic_bifurcation.png','-r300')