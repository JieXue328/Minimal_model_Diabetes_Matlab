function gluc_mm_mle
%GLUC_MM_MLE Maximum Likelihood Estimation of minimal model of glucose kinetics.
%History
%29-Jan-04, Natal van Riel, TU/e
%17-Jun-03, Natal van Riel, TU/e
%DATA
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
gluc_exp = tgi(:,2);
insul_exp = tgi(:,3);
%Fixed initial conditions
x0(2) = 0; %state variable denoting insulin action
%Fixed model parameters
Gb = gluc_exp(1); %[mg/dL] baseline glucose conc. in plasma
Ib = insul_exp(1); %[uU/mL] baseline insulin conc. in plasma
if 1 %normal subject
x0(1) = 279;%100; %[mg/dL] glucose conc. in plasma
Sg = 2.6e-2; %[1/min] glucose effectiveness
k3 = 0.025; %[1/min]
Si = 5.0e-4; %[mL/uU*min] insulin sensitivity
end
%REMARK 'if 1/0' CAN BE USED TO (DE)ACTIVATE THE CODE IN BETWEEN THE ifend
%STATEMENT
p = [Sg, Gb, k3, Si, Ib];
%Input
%insulin concentration in plasma [uU/mL]; assumed to be known at each simulation
%time sample from linear interpolation of its measured samples
tu = tgi(:,[1,3]);
t_insu = tu(:,1);
u = tu(:,2);
t_gluc = t_insu;
figure; plot(t_insu, u, 'o'); hold on
plot( [t_insu(1) t_insu(end)], [Ib Ib], '--k','Linewidth',1.5)
%baseline level
xlabel('t [min]'); ylabel('[\muU/mL]')
title('measured input signal (insulin conc. time course)')
tspan = [0:1:200]; %to verify interpolation of input signal
h = plot(tspan, interp1(tu(:,1),tu(:,2), tspan), '.r');
legend(h,'interpolated (resampled) signal')
%Simulation time vector:
tspan = t_insu;%tspan = union(t_insu, t_gluc);
disp(' Enter to continue with MLE'); disp(' ')
pause
if 1
%4 unknown model parameters:
p_init = [Sg %glucose effectiveness
k3 %
Si %insulin sensitivity
x0(1)]; %G0 initial glucose conc. in plasma
%3 known parameters:
p_fix = [Gb Ib x0(2)];
end
sigma_mu = 0;
sigma_nu = 0;
lb = 0*ones(size(p_init));%[0 0 0 0];
ub = [];%ones(size(p_init));%[];
options = optimset('Display','iter','TolFun', 1e-4,...%'iter' default:1e-4
'TolX',1e-5,... %default: 1e-4
'LevenbergMarquardt','on',... %default: on
'LargeScale','on'); %default: on
plt = 0;
%for i=1
%LSQNONLIN: objective function should return the model error
[p_est,J,RESIDUAL,exitflag,OUTPUT,LAMBDA,Jacobian] = lsqnonlin(@obj_fn,p_init,lb,ub,options,...
p_fix,gluc_exp,tspan,tu, sigma_nu,sigma_mu,plt); disp(' ')
%end
%Accuracy:
%lsqnonlin returns the jacobian as a sparse matrix
varp = resnorm*inv(Jacobian'*Jacobian)/length(tspan);
stdp = sqrt(diag(varp)); %The standard deviation is the square root
of the variance
%p = [Sg, Gb, k3, Si, Ib];
p = [p_est(1), p_fix(1), p_est(2), p_est(3), p_fix(2)];
%x0 = [G0, X0]
x0 = [p_est(4), p_fix(3)];
disp(' Parameters:')
disp([' Sg = ', num2str(p_est(1)), ' +/- ', num2str(stdp(1))])
disp([' Gb = ', num2str(p_fix(1))])
disp([' k3 = ', num2str(p_est(2)), ' +/- ', num2str(stdp(2))])
disp([' Si = ', num2str(p_est(3)), ' +/- ', num2str(stdp(3))])
disp([' Ib = ', num2str(p_fix(2))])
disp(' Initial conditions:')
disp([' G0 = ', num2str(p_est(4)), ' +/- ', num2str(stdp(4))])
disp([' X0 = ', num2str(p_fix(3))]); disp(' ')
plt = 1;
gluc = gluc_sim(tspan,x0,tu, p,sigma_nu,sigma_mu,plt);
%compare model output with measured data
figure(2); subplot(221)
h1 = plot(t_gluc,gluc_exp,'or', 'Linewidth',2);
%legend(h1, 'experimental test data')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function e = obj_fn(p_var, p_fix,data,tspan,tu, sigma_nu,sigma_mu,plt)
%4 unknown model parameters:
% p_var = [Sg k3 Si G0];
%3 known parameters:
% p_fix = [Gb Ib x0(2)];
%p = [Sg, Gb, k3, Si, Ib];
p = [p_var(1), p_fix(1), p_var(2), p_var(3), p_fix(2)];
%x0 = [G0, X0]
x0 = [p_var(4), p_fix(3)];
gluc = gluc_sim(tspan,x0,tu, p, sigma_nu,sigma_mu,0);
%LSQNONLIN: objective function should return the model error
e = gluc-data;
if plt==1 %fast update during estimation process
N=length(tspan);
figure(2)
plot(1:N,gluc,'b',1:N,data,'r'); drawnow
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gluc = gluc_sim(tspan,x0,tu, p, sigma_nu,sigma_mu,plt)
%SIM_INPUT_SIM Simulation of glucose minimal model.
%x0: initial conditions for glucose conc. and state variable denoting
% insulin action
%gluc: model output, glucose conc. in plasma [mg/dL]
if tspan(1) < tu(1,1)
error('no input defined at start of simulation (compare tspan(1))vs. tu(1,1)')
end
if 1 %Forward Euler integration
td = .2*min(diff(tspan)); %simulation time step is taken as 1/5 of
the
%minimal time interval between the samples of the experimental time
%vector
tspan_res = [];
for i = 1:length(tspan)-1
tspan_res = [tspan_res, tspan(i) :td: tspan(i+1)-td];
%resample while assuring that all time
%samples of tspan also occur
in tspan_res
end
tspan_res = [tspan_res, tspan(end)];
TD=diff(tspan_res); %vector with integration time steps (most equal to td)
u = interp1(tu(:,1),tu(:,2), tspan_res);
N = length(tspan_res);
x = x0;
for k = 2:N
dx(:,k) = gluc_ode([],x(k-1,:),u(k), p,sigma_nu); %column
x(k,:) = x(k-1,:) + dx(:,k)'*TD(k-1); %different
time samples in rows
end
for i = 1:length(tspan)
id(i) = find(tspan_res==tspan(i)); %select samples
corresponding to experiment
end
x = x(id,:);
elseif 0 %
ode_options = [];
[t,x] = ode45(@gluc_ode,tspan,x0,ode_options, tu, p,sigma_nu);
%[t,x] = ode15s(@gluc_ode,tspan,x0,ode_options, tu, p,sigma_nu);
end
%Output
gluc = x(:,1);
if plt==1
gluc_plt(tspan,x,tu,p,sigma_nu,sigma_mu)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dxout = gluc_ode(t,xin,tu, p,sigma_nu)
%GLUC_ODE ODE's of glucose minimal model
% xin: state values at previous time sample x(k-1) = [G(k-1); X(k-1)]
% dxout: new state derivatives dx(k) = [dG(k); dX(k)], column vector
% tu: u(k) (input signal at sample k) OR matrix tu
% p = [Sg, Gb, k3, Si, Ib];
%History
%17-Jun-03, Natal van Riel, TU/e
idG = 1; idX = 2; %glucose conc. [mg/mL] and state variable denoting
insulin action [uU/mL]
Sg = p(1); %[1/min] glucose effectiveness
Gb = p(2); %[mg/mL] baseline glucose conc.
k3 = p(3); %[1/min]
Si = p(4); %[mL/uU*min] insulin sensitivity
Ib = p(5); %[uU/mL] baseline insulin conc. in plasma
if length(tu)==1 %u(k) is provided as input to this function
u = tu;
else %calculate u(k) by interpolation of tu
u = interp1(tu(:,1),tu(:,2), t);
end
%[t u]
%ode's
dG = Sg*(Gb - xin(idG)) - xin(idX)*xin(idG);
dX = k3*( Si*(u-Ib) - xin(idX) );
dxout = [dG; dX];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gluc_plt(tspan,x,tu,p,sigma_nu,sigma_mu)
Gb = p(2);
Ib = p(5);
figure
subplot(221); h = plot(tspan,x(:,1), '-','Linewidth',2); hold on
plot( [tspan(1) tspan(end)], [Gb Gb], '--k','Linewidth',1.5)
%baseline level
ylabel('glucose level [mg/dL]')
legend(h, 'model simulation')
u = interp1(tu(:,1),tu(:,2), tspan); %reconstruct used input signal
subplot(222); plot(tspan,u, '--or','Linewidth',2); hold on
plot( [tspan(1) tspan(end)], [Ib Ib], '--k','Linewidth',1.5)
%baseline level
ylabel('insulin level [\muU/mL]'); xlabel('time [min]')
legend('interpolated measured test data')
%figure
subplot(223); plot(tspan, x(:,2), 'Linewidth',2)
xlabel('time [min]'); ylabel('interstitial insulin [1/min]')