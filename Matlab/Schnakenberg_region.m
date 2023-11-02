ccc
fs=15;
b=linspace(-1,2,1e3);
[x, y] = meshgrid(b);  % get 2-D mesh for x and y
conditions = (y< x/2) & (-(1/2)*x.^3+(1/2)*x < y);
cond = nan(length(b)); % Initialize
cond(conditions) = -1;

hold on
surf(x, y, cond)
shading interp
view(0,90)

p1=plot(b,b/2);
p2=plot(b,-(1/2)*b.^3+(1/2)*b);
plot(b,0*b,'k','linewidth',1)
xlabel('$\beta$')
ylabel('$\alpha$')
axis([0 2 -1 1])
L=legend([p1,p2],'$\beta$','$(-\beta^3+\beta)/2$','location','sw')
%set(L,'interpreter','latex')
set(gca,'fontsize',fs)

export_fig('../Pictures/Stable_alpha_region.png','-r300')
