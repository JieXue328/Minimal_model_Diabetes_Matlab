function gluc_mm_mle_trial
%GLUC_MM_MLE Maximum Likelihood Estimation of minimal model of glucose kinetics.

close all; clear all
%DATA
% data from Dr. Claudio Cobelli provided IGVTT data
%time (minutes) glucose level (mg/dl) insulin level (uU/ml)
tgi = [0	92.78	25.09	
2	217.99	357.90	
3	208.99	381.10	
4	208.99	338.40	
5	199.98	303.00	
8	179.62	130.10	
10	179.80	162.40	
18	161.78	83.37	
20	160.52	101.90	
23	159.80	90.63	
28	140.34	85.45	
32	139.08	78.43	
40	119.09	80.62	
60	96.03	46.37	
120	81.97	25.67	
180	78.37	14.68	
240	86.66	17.58];

gluc_exp = tgi(:,2);
insul_exp = tgi(:,3);
%Fixed initial conditions
x0(2) = 0; %state variable denoting insulin action
%Fixed model parameters
Gb = 82.51; %[mg/dL] baseline glucose conc. in plasma
Ib = 16.13; %[uU/mL] baseline insulin conc. in plasma

x0(1) = 279;%100; %[mg/dL] glucose conc. in plasma
Sg = 2.6e-2; %[1/min] glucose effectiveness
k3 = 0.025; %[1/min]
Si = 5.0e-4; %[mL/uU*min] insulin sensitivity

%STATEMENT
p = [Sg, Gb, k3, Si, Ib];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Input
%insulin concentration in plasma [uU/mL]; assumed to be known at each simulation
%time sample from linear interpolation of its measured samples
tspan = tgi(:,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

lb = 0*ones(size(p_init));%[0 0 0 0];
ub = [];%ones(size(p_init));%[];
options = optimset('Display','iter','TolFun', 1e-4,...%'iter' default:1e-4
'TolX',1e-5,... %default: 1e-4
'MaxFunEvals' , [],... %800 default:
'LargeScale','on'); %default: on
plt = 0;
%for i=1
%LSQNONLIN: objective function should return the model error
[p_est,resnorm,RESIDUAL,exitflag,OUTPUT,LAMBDA,Jacobian]= lsqnonlin(@err_fn,p_init,lb,ub,options,p_fix,tspan,tgi, plt); 
disp(' ')
%end
%Accuracy:
%lsqnonlin returns the jacobian as a sparse matrix
varp = resnorm*inv(Jacobian'*Jacobian)/length(tspan);
stdp = sqrt(diag(varp)); %The standard deviation is the square root of the variance
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
gluc = gluc_sim(tspan,x0,tgi, p,plt);
%compare model output with measured data
figure(1); subplot(221)
h1 = plot(tspan,gluc_exp,'or', 'Linewidth',2);
%legend(h1, 'experimental test data')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function e = err_fn(p_var, p_fix,tspan,tgi,plt)
p = [p_var(1), p_fix(1), p_var(2), p_var(3), p_fix(2)];
%x0 = [G0, X0]
x0 = [p_var(4), p_fix(3)];
gluc = gluc_sim(tspan,x0,tgi, p, 0);
%LSQNONLIN: objective function should return the model error
e = gluc-tgi(:,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gluc = gluc_sim(tspan,x0,tgi, p,plt)
%SIM_INPUT_SIM Simulation of glucose minimal model.
%x0: initial conditions for glucose conc. and state variable denoting
% insulin action
%gluc: model output, glucose conc. in plasma [mg/dL]

ode_options = [];
[t,x] = ode45(@gluc_ode,tspan,x0,ode_options, tgi, p);
%[t,x] = ode15s(@gluc_ode,tspan,x0,ode_options, tu, p,sigma_nu);

%Output
gluc = x(:,1);
if plt==1
gluc_plt(tspan,x,tgi,p)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dxout = gluc_ode(t,xin,tgi, p)
%GLUC_ODE ODE's of glucose minimal model
% xin: state values at previous time sample x(k-1) = [G(k-1); X(k-1)]
% dxout: new state derivatives dx(k) = [dG(k); dX(k)], column vector
% tu: u(k) (input signal at sample k) OR matrix tu
% p = [Sg, Gb, k3, Si, Ib];
%History
%17-Jun-03, Natal van Riel, TU/e

Sg = p(1); %[1/min] glucose effectiveness
Gb = p(2); %[mg/mL] baseline glucose conc.
k3 = p(3); %[1/min]
Si = p(4); %[mL/uU*min] insulin sensitivity
Ib = p(5); %[uU/mL] baseline insulin conc. in plasma

u = interp1(tgi(:,1),tgi(:,3), t);

%[t u]
%ode's
dG = Sg*(Gb - xin(1)) - xin(2)*xin(1);
dX = k3*( Si*(u-Ib) - xin(2) );
dxout = [dG; dX];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gluc_plt(tspan,x,tgi,p)
Gb = p(2);

Ib = p(5);
figure
subplot(221); h = plot(tspan,x(:,1), '-','Linewidth',2); hold on
plot( [tspan(1) tspan(end)], [Gb Gb], '--k','Linewidth',1.5)
%baseline level
ylabel('glucose level [mg/dL]')
legend(h, 'model simulation')

u = interp1(tgi(:,1),tgi(:,3), tspan); %reconstruct used input signal
subplot(222); plot(tspan,u, '--or','Linewidth',2); hold on
plot( [tspan(1) tspan(end)], [Ib Ib], '--k','Linewidth',1.5)
%baseline level
ylabel('insulin level [\muU/mL]'); xlabel('time [min]')
legend('interpolated measured test data')
%figure
subplot(223); plot(tspan, x(:,2), 'Linewidth',2)
xlabel('time [min]'); ylabel('interstitial insulin [1/min]')