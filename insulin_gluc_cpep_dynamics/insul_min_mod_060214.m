function insul_min_mod_060214
%INSUL_MIN_MOD Minimal model of insulin kinetics.

%Fixed initial conditions
%Fixed model parameters

Gb = 92;%118; %[mg/dL] baseline glucose conc. in plasma
Ib = 11;%10; %[uU/mL] baseline insulin conc. in plasma
y0 = 409.5;%100; %[uU/mL] insulin conc. in plasma at t0
k = 0.290; %[1/min] clearance rate of plasma insulin
gamma = 0.0055 % [1/min^2] meature of the secondary pancreatic response to glucose
p = [Gb, y0,k,gamma]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Input
%glucose concentration in plasma [mg/dL];
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

tg = tgi(:,[1,2]); % extract time-gluocose matrix
t_gluc = tg(:,1);  % time for gluocose matrix
g = tg(:,2);       % glucose matrix
t_insu = t_gluc;   % define t_insu
gluc_exp = tgi(:,2);    % denote glucose value matrix
insul_exp = tgi(:,3);   % denote insulin value matrix

figure; plot(t_gluc, g, 'o'); hold on                            % draw scatter plot for glucose-time experimental profile (in figure1)
plot( [t_gluc(1) t_gluc(end)], [Gb Gb], '--k','Linewidth',1.5)   % draw the baseline level for glucose(Gb)
xlabel('t [min]'); ylabel('glucose level [mg/dL]')
title('measured input signal (glucose conc. time course)')

tspan = 0:1:200; %to verify interpolation of input signal
h = plot(tspan, interp1(tg(:,1),tg(:,2), tspan), '.r');         % draw the interpolation for time-glucose experimental profile (in figure1)
legend(h,'interpolated (resampled) signal')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tspan = t_gluc; % reset tspan as actural experimental interval

ode_options = [];
[t,y] = ode45(@insu_ode,tspan,y0,ode_options,tg,p);
    
%Output
insu = y;  % get simulated insulin value matrix

figure
h = plot(tspan,insu,'-', 'Linewidth',2); hold on  % plot simulated time-insulin profile 
plot( [tspan(1) tspan(end)], [Ib Ib], '--k','Linewidth',1.5)    % plot the baseline level for insulin(Ib) (in figure2)
m = scatter (tspan,insul_exp,'o','r')
legend([h,m],'simulated time-insulin profile','experimental time-insulin profile')
ylabel('insulin level [uU/ml]')
title('INSULIN MINMAL MODEL')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dy = insu_ode(t,y,tg,p)
%INSU_ODE ODE's of insulin minimal model

Gb = p(1);
y0 = p(2);
k = p(3);
gamma = p(4);

g = interp1(tg(:,1),tg(:,2),t);

%ode's
dy = gamma*(g - Gb)*t - k * y; % xin(1) indicates G(t); xin(2) indicates X(t)


