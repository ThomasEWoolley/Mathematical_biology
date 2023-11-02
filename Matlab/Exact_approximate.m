ccc

D=10;
r=1;
L=linspace(0,10,1e2);
y1=sqrt(D/r)*acosh(0.5*cosh(sqrt(r/D)*L));
y2=L-sqrt(D/r)*log(2);

hold on
plot(L,real(y1),'b')
plot(L,y2,'r')
xlabel('$L$')
ylabel('$x_s$')
legend('Exact','Approximate','location','best')
export_fig('../Pictures/Exact_approximate.png','-r300')
