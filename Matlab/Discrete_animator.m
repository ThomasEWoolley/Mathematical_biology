ccc
fs=20;
a=0:0.001:4;
N=2000; x0=0.1;% Here I can write a comment
x=x0;
r=3.47;
Nits=40;
%% Bifurcation diagram
close all
figure('position',[0 0.1 1 1/2])
subplot(1,3,1)
hold on
for n=0:N
    x=a.*x.*(1-x);
    
    if n>.5*N
        plot(a,x,'b.','MarkerSize',3)
    end
    
end
axis ([0 4 0 1])
hold off
ylabel('$N^*$')
xlabel('$r$')
set(gca,'fontsize',fs)
hold on
plot([r r],[0 1],'r')
xticks(1:4)

%% Iterator
x=[];
x=x0;
for n=1:Nits
    x(n+1)=r.*x(n).*(1-x(n));
end

    xx=linspace(0,1);
for Iter=1:Nits
    %% Cobwebbing
    subplot(1,3,2)

    hold on
    plot(xx,xx,'r')
    plot(xx,r*xx.*(1-xx),'k')
    plot([x(1) x(1)],[0,x(2)],'b')
    for i=1:Iter
        plot([x(i) x(i)],[x(i) x(i+1)],'b')
        plot([x(i) x(i+1)],[x(i+1) x(i+1)],'b')
    end
    
    xlabel('$N_t$')
    ylabel('$N_{t+1}$')
    set(gca,'fontsize',fs)
    
    %% Plotter
    subplot(1,3,3)
    plot(0:length(x(1:Iter))-1,x(1:Iter),'b')
    ylabel('$N_t$')
    xlabel('$t$')
    axis([0 length(x)-1 0 1])
    set(gca,'fontsize',fs)
%     drawnow
    pause(0.1)
end