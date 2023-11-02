ccc
r=1;
x=linspace(0,6);
f=x.*exp(r*(1-x));
Nt(1)=3;
for i=1:100
    Nt(i+1)=Nt(i)*exp(r*(1-Nt(i)));
end
subplot(1,2,1)
hold on
plot(x,f,'r')
plot(x,x,'r')
% plot(Nt(1:end-1),Nt(2:end),'b')
axis equal
axis([0 6 0 2.5])
subplot(1,2,2)
plot(Nt)

%%
ccc
% fs=15;
figure('position',[0 0 1 .5])
r=1;
x=linspace(0,4,1e3);
f=x.*exp(r*(1-x));
hold on
plot(x,x,'r')
plot(x,f,'b')
formatter
export_fig('../Pictures/Ricker_r1.png','-r300')

figure('position',[0 0 1 .5])
r=2.1;
f=x.*exp(r*(1-x));
hold on
plot(x,x,'r')
plot(x,f,'b')
formatter
export_fig('../Pictures/Ricker_r21.png','-r300')

figure('position',[0 0 1 .5])
r=3;
f=x.*exp(r*(1-x));
hold on
plot(x,x,'r')
plot(x,f,'b')
formatter
export_fig('../Pictures/Ricker_r3.png','-r300')

function formatter
fs=20;
axis([0 4 0 2.5])
xlabel('$N_{t}$')
ylabel('$N_{t+1}$')
set(gca,'fontsize',fs)
end