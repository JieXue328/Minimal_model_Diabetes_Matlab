function ovallODE_060914

%Fixed initial conditions & model parameters

Gb = 92; %[mg/dL] baseline glucose conc. in plasma
Ib = 11;%[uU/mL] baseline insulin conc. in plasma
Cb = 54.57 %[uU/mL] baseline C-PEPTIDE conc. in plasma
x0(1) = 279;%[mg/dL] glucose conc. in plasma
x0(2) = 0; 
x0(3) = 409.5;%100; %[uU/mL] insulin conc. in plasma at t0
x0(4) = 159.827;
x0(5) = 1.28
Sg = 2.6e-2; %[1/min] glucose effectiveness
k3 = 0.025; %[1/min]
Si = 5.0e-4; %[mL/uU*min] insulin sensitivity
k = 0.290; %[1/min] clearance rate of plasma insulin
gamma = 0.0055 % [1/min^2] meature of the secondary pancreatic response to glucose
f = 0.818; 
kI = 0.122; 
kc = 0.028; 
p = [Gb,Ib,Sg, k3, Si,k,gamma,f, kI, kc]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% call ODE function
tspan = 0:5:180;
ode_options = [];
[t,x] = ode45(@all_ode,tspan,x0,ode_options,p);

gluc = x(:,1);  % get simulated glucose value matrix
X = x (:,2);    % get simulated interstitial insulin(X) value matrix
insul = x (:,3); % get simulated insulin value matrix
cpept = x (:,4);
r = x (:,5)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotting
figure
fig1 = plot(tspan,gluc,'r', 'Linewidth',2); hold on  % plot simulated time-glucose profile 
plot( [tspan(1) tspan(end)], [Gb Gb], '--k','Linewidth',1.5)    % plot the baseline level for glucose(Gb) (in figure3.1)
xlabel('time [min]'); ylabel('glucose level [mg/dL]')
title('GLUCOSE_TIME_PROFILE')

figure
fig2 = plot(tspan,insul,'g', 'Linewidth',2); hold on  % plot simulated time-insulin profile 
plot( [tspan(1) tspan(end)], [Ib Ib], '--k','Linewidth',1.5)    % plot the baseline level for insulin(Ib) (in figure2)
xlabel('time [min]'); ylabel('insulin level [uU/ml]')
title('INSULIN_TIME_PROFILE')

figure 
fig3 = plot(tspan, X,'b', 'Linewidth',2)                     % plot simulated time-interstitial insulin(X) profile
xlabel('time [min]'); ylabel('interstitial insulin [1/min]')
title('INTERSTITIAL INSULIN_TIME_PROFILE')

figure
fig4 = plot(tspan,c_pept,'p', 'Linewidth',2); hold on   
plot( [tspan(1) tspan(end)], [Cb Cb], '--k','Linewidth',1.5)    
xlabel('time [min]'); ylabel('C-peptide level [pmol/mL]')
title('C-PEPTIDE_TIME_PROFILE')

figure 
fig5 = plot(tspan, r,'c', 'Linewidth',2)                     
xlabel('time [min]'); ylabel('insulin secrete rate [1/min]')
title('INSULIN SCRETE RATE_TIME_PROFILE')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dxout = all_ode(t,xin,p)

Gb = p(1);
Ib = p(2);
Sg = p(3);
k3 = p(4);
Si = p(5);
k = p(6);
gamma = p(7);
f = p(8);
kI = p(9);
kc = p(10);

%ode's
dG = Sg*(Gb - xin(1)) - xin(1)*xin(2); % xin(1) indicates G(t); xin(2) indicates X(t); xin(3) indicates I(t)
dX = k3*( Si*(xin(3)-Ib) - xin(2) );
dI = gamma*(xin(1) - Gb)*t - k * xin(3); 
dI = f* xin(5) - kI*xin(3);          % xin(3) denotes I(t); xin(4) denotes C(t); xin(5) denotes r(t) 
dC = xin(5)  - kc* xin(4); 
dxout = [dG; dX;dI; dC; xin(5)];