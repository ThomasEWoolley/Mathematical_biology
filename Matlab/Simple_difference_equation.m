ccc
hold on
fs=15;
for r=[-1.5 -1 -0.5]% 0.5 1 1.5]
    n(1)=1;
    for i=1:5
        n(i+1)=r*n(i);
    end
    plot(0:5,n,'-o')
end
xlabel('$t$')
ylabel('$N_t$')
legend('$r=-1.5$','$r=-1$','$r=-0.5$','location','best')
set(gca,'fontsize',fs)
axis tight
export_fig('../Pictures/Simple_difference_equation_neg.png','-r300')

ccc
hold on
fs=15;
for r=[0.5 1 1.5]
    n(1)=1;
    for i=1:5
        n(i+1)=r*n(i);
    end
    plot(0:5,n,'-o')
end
xlabel('$t$')
ylabel('$N_t$')
legend('$r=0.5$','$r=1$','$r=1.5$','location','best')
set(gca,'fontsize',fs)
axis tight
export_fig('../Pictures/Simple_difference_equation_pos.png','-r300')
