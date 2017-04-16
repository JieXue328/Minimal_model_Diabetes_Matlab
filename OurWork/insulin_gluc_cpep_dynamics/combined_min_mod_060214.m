function combined_min_mod_060214

%Fixed initial conditions & model parameters

x0(2) = 0;

Gb = 92; %[mg/dL] baseline glucose conc. in plasma
Ib = 11;%[uU/mL] baseline insulin conc. in plasma
x0(1) = 279;%[mg/dL] glucose conc. in plasma
Sg = 2.6e-2; %[1/min] glucose effectiveness
k3 = 0.025; %[1/min]
Si = 5.0e-4; %[mL/uU*min] insulin sensitivity
x0(3) = 409.5;%100; %[uU/mL] insulin conc. in plasma at t0
k = 0.290; %[1/min] clearance rate of plasma insulin
gamma = 0.0055 % [1/min^2] meature of the secondary pancreatic response to glucose
p = [Gb,Ib,Sg, k3, Si,k,gamma]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Input
%time (minutes) glucose level (mg/dl) insulin level (uU/ml)
tgi = [ 0 92 11
    2 350 26
    4 287 130
    6 251 85
    8 240 51
    10 216 49
    12 211 45
    14 205 41
    16 196 35
    19 192 30
    22 172 30
    27 163 27
    32 142 30
    42 124 22
    52 105 15
    62 92 15
    72 84 11
    82 77 10
    92 82 8
    102 81 11
    122 82 7
    142 82 8
    162 85 8
    182 90 7];
tspan = tgi (:,1);      % denote time matix
gluc_exp = tgi(:,2);    % denote glucose value matrix
insul_exp = tgi(:,3);   % denote insulin value matrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ode_options = [];
[t,x] = ode45(@combined_ode,tspan,x0,ode_options,p);

%Output
gluc = x(:,1);  % get simulated glucose value matrix
X = x (:,2);    % get simulated interstitial insulin(X) value matrix
insul = x (:,3); % get simulated insulin value matrix

%plotting
figure
fig1 = plot(tspan,gluc,'r', 'Linewidth',2); hold on  % plot simulated time-glucose profile 
plot( [tspan(1) tspan(end)], [Gb Gb], '--k','Linewidth',1.5)    % plot the baseline level for glucose(Gb) (in figure3.1)
exp_gluc = scatter (tspan,gluc_exp,'o')
legend([fig1,exp_gluc],'simulated time-glucose profile','experimental time-glucose profile')
xlabel('time [min]'); ylabel('glucose level [mg/dL]')
title('GLUCOSE_TIME_PROFILE')

figure
fig2 = plot(tspan,insul,'g', 'Linewidth',2); hold on  % plot simulated time-insulin profile 
plot( [tspan(1) tspan(end)], [Ib Ib], '--k','Linewidth',1.5)    % plot the baseline level for insulin(Ib) (in figure2)
exp_insul = scatter (tspan,insul_exp,'o')
legend([fig2,exp_insul],'simulated time-insulin profile','experimental time-insulin profile')
xlabel('time [min]'); ylabel('insulin level [uU/ml]')
title('INSULIN_TIME_PROFILE')

figure 
fig3 = plot(tspan, X,'b', 'Linewidth',2)                     % plot simulated time-interstitial insulin(X) profile
legend(fig3,'simulated time-interstitial insulin profile')
xlabel('time [min]'); ylabel('interstitial insulin [1/min]')
title('INTERSTITIAL INSULIN_TIME_PROFILE')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dxout = combined_ode(t,xin,p)
%COMBINED_ODE ODE's of glucose/insulin minimal model

Gb = p(1);
Ib = p(2);
Sg = p(3);
k3 = p(4);
Si = p(5);
k = p(6);
gamma = p(7);

%ode's
dG = Sg*(Gb - xin(1)) - xin(1)*xin(2); % xin(1) indicates G(t); xin(2) indicates X(t); xin(3) indicates I(t)
dX = k3*( Si*(xin(3)-Ib) - xin(2) );
dI = gamma*(xin(1) - Gb)*t - k * xin(3); 
dxout = [dG; dX;dI];


