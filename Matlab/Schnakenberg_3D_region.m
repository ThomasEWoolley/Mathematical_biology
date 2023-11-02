ccc
figure('position',[0 0 1/4 .4])
fs=15;
n=1e2;
b1=linspace(0,2,10*n);
b2=linspace(0,50,n);
[x,y,z] = meshgrid(b1,b1,b2);  % get 2-D mesh for x and y
conditions = (y< x/2) & (-(1/2)*x.^3+(1/2)*x < y) & ((-2*y+3*x+2*sqrt(-2*x.*y+2*x.^2)).*x.^3./(4*y.^2-4*y.*x+x.^2)<z);


hold on
p=patch(isosurface(x,y,z,conditions,.5));
p.FaceColor = 'cyan';
p.EdgeColor = 'none';
view([30 20]); 
camlight(30,20)
xlabel('$\beta/2$')
ylabel('$\alpha$')
zlabel('$D_z$')
p1=plot(b1,b1/2);
p2=plot(b1,-(1/2)*b1.^3+(1/2)*b1);
axis([0 2 0 1 0 50])
L=legend([p1,p2,p],'$\beta/2$','$(-\beta^3+\beta)/2$','$D_v>D_{+v}$','location','best');
set(gca,'fontsize',fs)
export_fig('../Pictures/Stable_patterning_region.png','-r300')
