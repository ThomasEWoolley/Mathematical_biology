function Logistic(varargin)
% Simple ODE solver for the logistic equation,
% du/dt=r*u*(c-u), u(0)=u0.
% If the function is called without parameters then default parameters,
% c=10, r=0.1, u0=1
% are used.
% Alternatively, these three parameters can be passed to the function
% in the form
% Logistic(c,r,u0).
% (c) Thomas E. Woolley, 2017.


% Clears screen.
clc
close all

cont=1; % Checking variable.

% The if loop checks how many variables are specified.
if isempty(varargin)
     c=10;
     r=0.1;
     u0=0.1;
elseif length(varargin)~=3
    disp('This function must have either 3, or no, parameters specified')
    cont=0;
else    
     c=varargin{1};
     r=varargin{2};
     u0=varargin{3};
end

if cont==1
tspan = [0 10]; % Length of time the equation is solved for.
[t,u] = ode23(@(t,u) ODE(c,r,t,u), tspan, u0); % Matlab ODE solver.
plot(t,u,'linewidth',3) % Plot solution.
xlabel('Time')
ylabel('Population')
shg
end


function dydt = ODE(c,r,t,u)
% ODE equation to solve.
dydt = r*u*(c-u); % The Logistic equation.