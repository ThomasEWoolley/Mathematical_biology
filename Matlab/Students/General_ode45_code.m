%% Ensure we start from a blank slate
clear all
close all
clc

%% Initialise variables
fs=15; % Set fontsize

Tend=0.3; % Set final time
N=100; % Set number of steps
t=linspace(0,Tend,N); % Set number of time points at which the solution is to be evaluated
IC=[200;
    190]; %Initial conditions as a column vector

%% Solve ODE

% Each component is a row in the matrix y.
[t,y] = ode45(@ODE_function,t,IC);

%% Plot Solution

%%% (time,solution plots)
hold on % Allows us to plot multiple graphs at once
plot(t,y(:,1)) % Plots the first solution
plot(t,y(:,2)) % Plots the second solution
legend('u','v') % Adds a legend and defines the lines for clarity

% Label axes for clarity
xlabel('t') 
ylabel('Solutions')

set(gca,'fontsize',fs) % Changes the axis fontsize

%%% Phasespace plots
figure % Starts a new figure without overwriting the previous one
plot(y(:,1),y(:,2))
xlabel('u') 
ylabel('v')
set(gca,'fontsize',fs)




%% ODE description
function dydt=ODE_function(t,y)
r=1; %Interaction parameter

% Defining the ODEs as a column vector. Each row is one of the ODEs.
dydt=[-r*y(1)*y(2)
      -r*y(1)*y(2)]; 

end