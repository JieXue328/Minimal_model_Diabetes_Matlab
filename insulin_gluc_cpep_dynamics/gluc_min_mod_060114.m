function gluc_min_mod_060114
%GLUC_MIN_MOD Minimal model of glucose kinetics.

%Fixed initial conditions
%Fixed model parameters

x0(2) = 0;

Gb = 92;%118; %[mg/dL] baseline glucose conc. in plasma
Ib = 11;%10; %[uU/mL] baseline insulin conc. in plasma
x0(1) = 279;%100; %[mg/dL] glucose conc. in plasma
Sg = 2.6e-2; %[1/min] glucose effectiveness
k3 = 0.025; %[1/min]
Si = 5.0e-4; %[mL/uU*min] insulin sensitivity
p = [Sg, Gb, k3, Si, Ib]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Input
%insulin concentration in plasma [uU/mL];
%Data from Pacini & Bergman (1986) Computer Methods and Programs in Biomed. 23: 113-122.
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

tu = tgi(:,[1,3]); % extract time-insulin matrix
t_insu = tu(:,1);  % time for insulin matrix
u = tu(:,2);       % insulin matrix
t_gluc = t_insu;   % define t_gluc
gluc_exp = tgi(:,2);    % denote glucose value matrix
insul_exp = tgi(:,3);   % denote insulin value matrix

figure; subplot(221); plot(t_insu,insul_exp,'o', 'Linewidth',2); % draw scatter plot for insulin-time experimental profile (in figure1.1)
hold on
plot( [t_insu(1) t_insu(end)], [Ib Ib], '--k','Linewidth',1.5)   % draw the baseline level for insulin(Ib)
ylabel('insulin level [\muU/mL]'); xlabel('time [min]')

subplot(222); plot(t_gluc,gluc_exp, 'o','Linewidth',2); hold on  % draw scatter plot for glocose-time experimental profile  (in figure1.2)
plot( [t_gluc(1) t_gluc(end)], [Gb Gb], '--k','Linewidth',1.5)   % draw the baseline level for glucose(Gb)
ylabel('glucose level [mg/dL]'); xlabel('time [min]')

figure; plot(t_insu, u, 'o'); hold on                            % draw scatter plot for insulin-time experimental profile (in figure2)
plot( [t_insu(1) t_insu(end)], [Ib Ib], '--k','Linewidth',1.5)   % draw the baseline level for insulin(Ib)
xlabel('t [min]'); ylabel('[\muU/mL]')
title('measured input signal (insulin conc. time course)')

tspan = 0:1:200; %to verify interpolation of input signal
h = plot(tspan, interp1(tu(:,1),tu(:,2), tspan), '.r');         % draw the interpolation for time-insulin experimental profile (in figure2)
legend(h,'interpolated (resampled) signal')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tspan = t_insu; % reset tspan as actural experimental interval
u = interp1(tu(:,1),tu(:,2), tspan);  %reconstruct used input signal (define I(t))

ode_options = [];
[t,x] = ode45(@gluc_ode,tspan,x0,ode_options,tu,p);
    
%Output
gluc = x(:,1);  % get simulated glucose value matrix
X = x (:,2);    % get simulated interstitial insulin(X) value matrix
figure
subplot(221); h = plot(tspan,gluc,'-', 'Linewidth',2); hold on  % plot simulated time-glucose profile 
plot( [tspan(1) tspan(end)], [Gb Gb], '--k','Linewidth',1.5)    % plot the baseline level for glucose(Gb) (in figure3.1)
ylabel('glucose level [mg/dL]')
title('normal FSIGT')

subplot(222); plot(tspan,u,'--or', 'Linewidth',2); hold on      % plot the interpolation for time-insulin experimental profile (in figure3.1)
plot( [tspan(1) tspan(end)], [Ib Ib], '--k','Linewidth',1.5)    % plot the baseline level for insulin(Ib)
ylabel('insulin level [\muU/mL]');
xlabel('time [min]')
legend('interpolated measured test data')

subplot(223); plot(tspan, X, 'Linewidth',2)                     % plot simulated time-interstitial insulin(X) profile 
xlabel('time [min]'); ylabel('interstitial insulin [1/min]')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dxout = gluc_ode(t,xin,tu,p)
%GLUC_ODE ODE's of glucose minimal model

Sg = p(1);
Gb = p(2);
k3 = p(3);
Si = p(4);
Ib = p(5);

u = interp1(tu(:,1),tu(:,2),t);

%ode's
dG = Sg*(Gb - xin(1)) - xin(1)*xin(2); % xin(1) indicates G(t); xin(2) indicates X(t)
dX = k3*( Si*(u-Ib) - xin(2) );
dxout = [dG; dX];

