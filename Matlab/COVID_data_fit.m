ccc
fitVirusCV19v3('UK')


function res = fitVirusCV19v3(country,varargin)

global C dC
global w1 w2  % optimization weights
%minimal check
narginchk(1,inf)
nargoutchk(0,1)
% default values
Nmax = 2e6;   % max. population size
plt = true;   % plot
w1 = [];      % value weight
w2 = [];      % derivative weight
nsp = 2;

result=readtable('UK_confirmed.csv');
deathresult=readtable('UK_deaths.csv');
writetable(result,'result.txt','WriteVariableNames',false);
writetable(deathresult,'deathresult.txt','WriteVariableNames',false);
opts = detectImportOptions('result.txt', "TextType","string");
opts1 = detectImportOptions('deathresult.txt', "TextType","string");
first_day = datetime(2020,1,22)
day_add = size(result);
last_day = first_day+days(day_add(2)-5)
time = first_day:last_day;
C = cell(1,day_add(2));
for i = 1:day_add(2)
    if i == 1
        C(i) = {'Province_State'};
    elseif i == 2
        C(i) = {'Country_Region'};
    elseif i == 3
        C(i) = {'Lat'};
    elseif i == 4
        C(i) = {'Long'};
    else
        formatOut = 'xmm_dd_yy';
        C(i) = {sprintf('%s',datestr(datenum(time(i-4)),formatOut))};
    end
end
times_conf = readtable('result.txt',opts);
times_conf1 = readtable('deathresult.txt',opts1);
matlab.lang.makeValidName(C);
times_conf.Properties.VariableNames = C;
times_conf1.Properties.VariableNames = C;
times_conf.("Country_Region")(times_conf.("Country_Region") == "United Kingdom") = "UK";
times_conf1.("Country_Region")(times_conf1.("Country_Region") == "United Kingdom") = "UK";
vars = times_conf.Properties.VariableNames;
vars1 = times_conf1.Properties.VariableNames;
times_conf_country = groupsummary(times_conf,"Country_Region",{'sum'},vars(3:end));
times_conf_country1 = groupsummary(times_conf1,"Country_Region",{'sum'},vars1(3:end));
vars = times_conf_country.Properties.VariableNames;
vars = regexprep(vars,"^(sum_)(?=L(a|o))","remove_");
vars = erase(vars,{'sum_'});
times_conf_country.Properties.VariableNames = vars;
vars1 = times_conf_country1.Properties.VariableNames;
vars1 = regexprep(vars1,"^(sum_)(?=L(a|o))","remove_");
vars1 = erase(vars1,{'sum_'});
times_conf_country1.Properties.VariableNames = vars1;
infectedtable = removevars(times_conf_country,[{'GroupCount'},vars(contains(vars,"remove_"))]);
countrytable = infectedtable(strcmp(infectedtable.("Country_Region"),country), :);
deathtable = removevars(times_conf_country1,[{'GroupCount'},vars1(contains(vars1,"remove_"))]);
countrytable1 = deathtable(strcmp(deathtable.("Country_Region"),country), :);
countrytable = countrytable(:,2:end);
countrytable1 = countrytable1(:,2:end);
cols1 = size(countrytable);
cols2 = size(countrytable1);
Countrytotaldead = zeros(1,cols2(2));
Countrytotalinfected = zeros(1,cols1(2));
for i = 1:cols1(2)
    Countrytotalinfected(i) = table2array(countrytable(1,i));
    Countrytotaldead(i) = table2array(countrytable1(1,i));
end
C = Countrytotalinfected;
startidx = zeros(length(C));
for i = 1:length(C)
    if i == length(C)
        break
    end
    if C(i+1) > 1.5*C(i) && C(i) > 10
        startidx(i) = i;
    end
end
startidx = find(startidx~=0, 1, 'first');
first_day=datenum('2020/01/22'); % spread start date - do not change
date0 = first_day+startidx;
C = C(startidx:end);
% initial guess
n0 = 1;
C = C(n0:end);
b0 = iniGuess(C);

% ... logistic curve parameters
K0 = b0(1);
r  = b0(2);
A  = b0(3);
C0 = K0/(A + 1);
% ... initial guess
I0 = C0;
N = 2*K0;
gamma = 2*r;
beta  = 1.5*gamma;
% main calculation =======================================================%
% set infection rate and time intervals
dC  = diff(C);
dC(dC<0) = 0;  % correct
nday = length(C);
tt   = 0:nday-1;  % time span
% initial estimate
b0 = [beta gamma N I0]';
% calculate parameters
if ~isempty(w1) && ~isempty(w2)
    % weigts are set by user
    [b,~,~] = parest(b0);
else
    % automatic selection of weigths
    for i = 1:3
        switch i
            case 1
                w1 = 1;
                w2 = 1;
            case 2
                w1 = 1;
                w2 = 0;
            case 3
                w1 = 0;
                w2 = 1;
        end
        [b,~,~] = parest(b0);
        if all(b > 0) && b(3) <= Nmax
            break
        end
    end
end

% unpack results
beta  = b(1);
gamma = b(2);
N     = b(3);
I0    = b(4);
% postprocessing ======================================================== %
%... final value
Clim = calcClim(b);
%... value at inflection point
Cm   = calcCm(b);
% basic reproduction number
R0 = beta/gamma*(1 - I0/N);
%... parameters of logistic model approximation
r = beta - gamma;
K = 2*(beta - gamma)/(2*beta - gamma)*N;
t2 = log(2)/r;
%... tangent slope in inflection point
k = (N - Cm)*(beta*Cm/N + gamma*log((N - Cm)/(N - I0)));
%... acceleration time
tau1  = Cm/k;
%... deceleration time
tau2  = (Clim - Cm)/k;
%... total duration of accelerated phase
tau = tau1 + tau2;
%... inflection time
tm = calcTm(b,Cm);
tm = real(tm);
%... datums
tp1 = (tm - tau1) + date0;  % begin acceleration
tp2 = (tm) + date0;         % turning point
tp3 = (tm + tau2) + date0;  % end deceleration
tp4 = (tm + tau2) + tau + date0; % enter final phase
%... dense forcast curve
dt = 0.1;
tspan = 0:dt:2.5*tm;
warning('on')
[t,Ca] = ode45(@(t,y) odeFun(t,y,b), tspan, I0);
warning('off')
Ca = real(Ca);
%... calculate forcasting curve at data points
tspan = 0:nday; % one day more
warning('on')
[~,Ce] = ode45(@(t,y) odeFun(t,y,b), tspan, I0);
warning('off')
Ce = real(Ce);
Cnxt = Ce(end);   % one day forcast
if Cnxt < C(end)
    % model fails. Cnxt can not be less than current actual.
    Cnxt = NaN;
end
Ce =Ce(1:end-1);  % delete last
%... calculate statistics
[R2,~,RMSE,~,~] = calcR2(C',Ce(1:nday));
if R2 < 0.9
    fprintf('***Warning: R2 = %g\n',R2)
end


% plot results ===========================================================%

figure('position',[0 0 0.75 0.9])
%...set scale
if max(Ca) > 1000
    sf = 1000;
else
    sf = 1;
end
ttt  = 0:nday-2;
% plot total cases ---------------------
subplot(nsp,1,2)
hold on
%... plot forcast curve
plot(t + date0,Ca/sf,'LineWidth',2)
%... plot +/-SDE
nsd = 3;
% h = plot(t + date0,(Ca + nsd*RMSE)/sf,'r','LineWidth',1);
% h.Annotation.LegendInformation.IconDisplayStyle = 'off';
% Ct = (Ca - nsd*RMSE)/sf;
% Ct(Ct<0) = 0;
% h = plot(t + date0,Ct,'r','LineWidth',1);
% h.Annotation.LegendInformation.IconDisplayStyle = 'off';
%...get plot limits
ylm = get(gca,'Ylim');  % get y-axes limits
xlm = get(gca,'Xlim');  % get x-axes limits
www = xlm(2);
hhh = ylm(2);
%... plot cases limits
h = plot([0,t(end)] + date0,[Clim,Clim]/sf,'g--','LineWidth',1);
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
%... add data points
scatter(tt + date0, C/sf,50,'k','filled')
% h = scatter(tt + date0, C/sf,30,'w','filled');
% h.Annotation.LegendInformation.IconDisplayStyle = 'off';
%... limits
xlim([t(1),t(end)]+date0);
%... what kind of thicks?
datetick('x',19,'keeplimits')
%... label axes
xlabel('Date')
if sf == 1000
    ylabel('Infected + Recovered ($\times$ 1000 cases)','interpreter','latex')
else
    ylabel('Infected (cases)')
end
    fprintf('  Contact rate (beta)           %g (1/(person x day))\n',(beta/N))
    fprintf('  Removal rate (gamma)          %g (1/day)\n',round(gamma,3))
    fprintf('  Population size (N)           %g\n',fix(N))
    fprintf('  Initial number of cases (I0)  %g\n',fix(I0))
    fprintf('Contact number            (R )  %g\n',round(beta/N*1/gamma*(N-I0),3));
    fprintf('  Final number of cases         %g\n',fix(Clim))

%... add legend
legend('Prediction','Actual','Location','best')
set(gca,'fontsize',15)

% plot infection rate -------------
subplot(nsp,1,1)
hold on
%... plot data
bar( date0 + ttt,dC)
%... plot forcasting curve
plot(t(1:end-1) + date0,diff(Ca)/dt,'LineWidth',2)
%... limits
xlim([t(1),t(end)]+date0);
%... what kind of thicks?
datetick('x',19,'keeplimits')
% ... add labels
ylabel('New cases/day','interpreter','latex')
xlabel('Date')
%... add legend
legend('Actual','Prediction','Location','best')
set(gca,'fontsize',15)
export_fig('../Pictures/Covid.png','-r300')
end
function [b,fmin,flag] = parest(b0)
%PAREST Parameter estimation
%
%   This function use MATLAB's fminsearch
%
warning('on')
options = optimset('Display','off'); %,'MaxIter',maxiter,...
%  'MaxFunEvals',maxfun);
[b, fmin,flag] = fminsearch(@fun, b0, options);
warning('off')
%     fprintf('Exit condition: %g\n',flag)
%     fprintf('Smallest value of the error: %g\n',fmin);
end
function f = fun( par)
%FUN Optimization function
global C dC
global w1 w2
% upack parameter
I0 = par(4);
% set time span
tspan = 0:length(C)-1;
% solve ODE
try
    warning('off')
    [tsol,Csol] = ode45(@(t,y) odeFun(t,y,par), tspan, I0);
    warning('on')
catch
    f = NaN;
    warning('on')
    return
end
% check if calculation time equals sample time
if length(tsol) ~= length(tspan)
    f = NaN;
    return
end
% calculate optimization function
c1 = w1/(w1 + w2);
c2 = w2/(w1 + w2);
f1 = 0;
f2 = 0;
if c2 > 0
    f2 = norm((dC' - diff(Csol)));
end
if c1 > 0
    f1 = norm((C' - Csol));
end
f =  c1*f1  +  c2*f2;
end
function res = calcClim(par)
%CALCCLIM Calculate number of recoverd individuals after t=inf
beta  = par(1);
gamma = par(2);
N     = par(3);
I0    = par(4);
res = calcEndPoint(beta,gamma,I0/N)*N;
end
function res = calcCm(par)
%CALCCM Calculate number of cases at inflection point
beta  = par(1);
gamma = par(2);
N     = par(3);
I0    = par(4);
res = calcInflectionPoint(beta,gamma,I0/N)*N;
end
function res = calcTm(par,Cm)
%CALCTM Calculate peak time
beta  = par(1);
gamma = par(2);
N     = par(3);
c0    = par(4)/N;
warning('off')
res = integral(@fun,c0,Cm/N);
warning('on')
    function t = fun(c)
        tt = (1 - c).*(beta*c + gamma*log((1 - c)/(1 - c0)));
        t = 1./tt;
    end
end
function ce = calcEndPoint(beta,gamma,c0)
%CALCENDPOINT Calculate end density
ce = 1 + gamma/beta*...
    lambertw(-beta*(1 - c0)*exp(-beta/gamma)/gamma);
end
function cm = calcInflectionPoint(beta,gamma,c0)
%CALCINFLECTIONPOINT Calculate inflection point for density curve
cm = 1 + (gamma/2/beta)*...
    lambertw(-1, -2*beta*(1 - c0)*exp(-(1 + beta/gamma))/gamma);
end
function [R2, AdjR2, RMSE, Fval,pval] = calcR2(y,ye)
%CALCR2 Calculate the coefficient of determination
%
% Input:
%   y  -- actual values
%   ye -- estimated values
%
% Output:
%   R2 -- the coefficient of determination
%   AdjR2 -- adjusted R2
%
% References:
%   https://en.wikipedia.org/wiki/Coefficient_of_determination
%
n = length(y);  % number of data points
p = 4;          % number of explanatory terms in a model
ybar  = sum(y)/n;
SStot = sum((y - ybar).^2);
SSres = sum((y - ye).^2);    %
R2    = 1 - SSres/SStot;
% calculate adjusted R2
if nargout > 1
    AdjR2 = 1 - (1 - R2)*(n - 1)/(n - p - 1);
end
if nargout > 2
    % http://facweb.cs.depaul.edu/sjost/csc423/documents/f-test-reg.htm
    SSM =  sum((ye - ybar).^2); % sum of squares for regression
    SSE = SSres;                % sum of squares for residuals
    SST = SStot;                % sample variance x (n-1)
    dfm = p - 1;            % Corrected Degrees of Freedom for Model
    dfe = n - p;            %Degrees of Freedom for Error
    dft = n - 1;            %Corrected Degrees of Freedom Total
    MSM = SSM/dfm;          % Mean of Squares for Model
    MSE = SSE/dfe;          % Mean of Squares for Error (variance of the residuals)
    %MST = SST/DFT;         % Mean of Squares Total (sample variance)
    RMSE = sqrt(MSE); %sqrt(SStot/(n - p));  % standard error of estimate
    % calculate F statistics
    Fval = MSM/MSE;
    pval = fcdf(1/max(0,Fval),dfe,dfm); %????
end
end
function [b0] = iniGuess(C)
%INIGUESS Initial guess for logistic regression
% calculate initial K, r, A using data from three equidistant points
%
% Input:
%   C -- data
%
% Output:
%   b0 -- initial guess = [K r A]' or [] if calculation fails
b0 = [];
n = length(C);

if n <= 5
    fprintf('***Warning: not enough data.\n')
    return
end
nmax = n - 5;

for i = 1:nmax
    % calculate time interval for equidistant points: k-2*m, k-m, k
    if mod(n-i+1,2) == 0
        k1 = i;
        k3  = n-1;
    else
        k1 = i;
        k3 = n;
    end
    k2 = (k1 + k3)/2;
    m = k2 - k1 -1;
    if k1 <1 || k2 < 1 || k3 < 1 || m < 1
        break
    end
    
    if isnan(C(k1)) || isnan(C(k2)) || isnan(C(k3))
        continue
    end
    % calculate K, r, A ...
    %.. calculate K
    q = C(k2)^2 - C(k3)*C(k1);
    if q <= 0
        %      fprintf('***Warning: iniGuess q = %g  k1 = %d k2= %d k3 = %d \n',...
        %          q, k1, k2, k3)
        continue
    end
    p = C(k1)*C(k2) - 2*C(k1)*C(k3) + C(k2)*C(k3);
    if p <= 0
        %   fprintf('***Warning: iniGuess p = %g\n',p)
        continue
    end
    K = C(k2)*p/q;
    % ... calculate r
    r = log(C(k3)*(C(k2) - C(k1))/C(k1)/(C(k3) - C(k2)))/m;
    if r < 0
        %  fprintf('***Warning: iniGuess r = %g\n',r)
        continue
    end
    %... calculate A
    A = (C(k3) - C(k2))*(C(k2) - C(k1))/q*...
        (C(k3)*(C(k2) - C(k1))/C(k1)/(C(k3) - C(k2)))^((k3-m)/m);
    if A <= 0
        %   fprintf('***Warning: iniGuess A = %g\n',r)
        continue
    end
    % this is initial guess
    b0 = [K r A]';
    break
end
end
function dCdt = odeFun(~,C,par)
%ODEFUN SIR model
% unpack parameters
beta  = par(1);
gamma = par(2);
N     = par(3);
I0    = par(4);
% set temp. vars
c0    = I0/N;
c     = C/N;
% setup equation
dCdt = N*(1 - c)*(beta*c + gamma*log((1 - c)/(1 - c0)));
end
