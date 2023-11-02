%% Ensure we start from a blank slate
clear all
close all
clc

%% Initialise variables
fs=15; % Set fontsize

b=4;
E=0.8;

Np=(b+sqrt(b^2-4*(1+E)^2))/(2*(1+E)); %N+ from question 2.
Nn=(b-sqrt(b^2-4*(1+E)^2))/(2*(1+E)); %N- from question 2.
N1=(b-sqrt(b^2-4*E*(1+Nn^2)*(E*(1+Nn^2)-b*Nn)))/(2*E*(1+Nn^2)); %N1 from question 2.
N2=(b+sqrt(b^2-4*E*(1+Nn^2)*(E*(1+Nn^2)-b*Nn)))/(2*E*(1+Nn^2)); %N2 from question 2.

%% Set up the basic plotting space
N=linspace(-1,8);

hold on
plot(N,N,'k') % Plot N(t)=N(t+1)
plot(N,b*N.^2./(1+N.^2)-E*N,'b') % Plot N(t+1)=f(N(t))
plot(N,Nn*ones(1,length(N)),'--k','linewidth',1) % Plot N(t+1)=N-
plot([N1 N1],[0 Nn],'--k','linewidth',1) % Plot N(t)=N1
plot([Nn Nn],[0 Nn],'--k','linewidth',1) % Plot N(t)=N-
plot([N2 N2],[0 Nn],'--k','linewidth',1) % Plot N(t)=N2
axis([0 5 0 2])
xlabel('$N_t$')
ylabel('$N_{t+1}$')

%% Calculate the cobweb diagram
% The Cobweb(x1,x2,x3,x4) function has 4 arguments x1-x4.
% x1 is the initial point from which the cobweb starts.
% x2 and x3 are the parameters b and E, respectively.
% x4 is the colour of the cobweb.
Cobweb(3,b,E,'b')
Cobweb(0.5,b,E,'r')
Cobweb(4.5,b,E,'r')
set(gca,'fontsize',fs) % Set fontsize.

function Cobweb(N0,b,E,c)
% The following for loop computes the first 100 iteration values.
Nt(1)=N0;
for i=1:100
    Nt(i+1)=b*Nt(i)^2/(1+Nt(i)^2)-E*Nt(i);
end


% The following code plots the cobweb.
plot([Nt(1) Nt(1)],[0 Nt(2)],c,'linewidth',1)
for i=1:100-2
    plot([Nt(i) Nt(i+1)],[Nt(i+1) Nt(i+1)],c,'linewidth',1)
    plot([Nt(i+1) Nt(i+1)],[Nt(i+1) Nt(i+2)],c,'linewidth',1)
end
end