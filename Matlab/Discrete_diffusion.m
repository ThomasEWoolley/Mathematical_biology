%% Ensure we sart from a blank slate
clear all
close all
clc

%% Initialise variables
fs=15; % Set fontsize

dx=0.1; % Set box size
T_steps=1000; % Set number of time steps over which the solution is evaluated
x=0:dx:1; % Create the spatial discretisation
N=length(x); % Defines how many ODEs there will be
t=linspace(0,1,T_steps); % Create time series
ICs=zeros(N,1); % Intialise initial condition to be a zero vector
ICs(1)=2; % Put a value of 2 in the first box as the initial condition

%% Solve the N coupled ODEs
[t,y]=ode45(@(t,y)ODE(t,y,N,dx),t,ICs);


%% Plot space-time graph
figure('units','normalized','position',[0 0 1 1/3])
subplot(1,2,1)
pcolor(x,t,y)
shading interp
xlabel('Space')
ylabel('Time')
colorbar
set(gca,'fontsize',fs)

%% Plot spatial profile over a number of time points
subplot(1,2,2)
plot(x,y(1,:))
hold on
plot(x,y(floor(T_steps/10),:))
plot(x,y(T_steps,:))
legend('$t=0$',['$t=$ ',num2str(t(floor(T_steps/10)))],['$t=$ ',num2str(t(floor(T_steps)))])
xlabel('Space')
ylabel('$u$')
set(gca,'fontsize',fs)






function dydt=ODE(t,y,N,dx)
D=0.1; % Diffusion rate

% First and last box are set with Neumann boundary conditions
dydt(1,1)=D/dx^2*(-y(1)+y(2));
dydt(N,1)=D/dx^2*(y(N-1)-y(N));

% Fill in all other vector values
for i=2:N-1
    dydt(i,1)=D/dx^2*(y(i-1)-2*y(i)+y(i+1));
end

end