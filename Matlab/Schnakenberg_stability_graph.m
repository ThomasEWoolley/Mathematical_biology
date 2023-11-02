ccc
beta=linspace(0,3);
fp=(1-beta.^2+sqrt((1-beta.^2).^2-4*beta.^2))./2;
fn=(1-beta.^2-sqrt((1-beta.^2).^2-4*beta.^2))./2;

hold on
plot(beta,real(fp),'r')
plot(beta,imag(fp),'r--')
plot(beta,real(fn),'b')
plot(beta,imag(fn),'b--')

axis([0 3 -5 3])
xlabel('$\beta$','interpreter','latex')
ylabel('Eigenvalues','interpreter','latex')
l=legend('$\textrm{Re}(\lambda_+)$','$\textrm{Im}(\lambda_+)$','$\textrm{Re}(\lambda_-)$','$\textrm{Im}(\lambda_-)$','location','best');
set(l,'interpreter','latex')
set(gca,'fontsize',15)
export_fig('..\Pictures\Schnakenberg_stability.png','-r300')