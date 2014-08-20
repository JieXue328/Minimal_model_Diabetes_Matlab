function Cpep_oral_mm_mle_fixed_081414

tgc = [0	5.39	0.62
10	5.93	0.84
20	6.95	1.27
30	8.05	1.68
45	9.04	2.11
60	9.44	2.27
75	9.27	2.41
90	8.84	2.38
105	8.25	2.39
120	7.58	2.28
135	7.05	2.19
150	7.03	2.09
165	6.54	1.93
180	6.01	1.7
195	5.45	1.45
210	5.04	1.26
225	4.78	1.17
240	4.48	0.93
255	4.34	0.76
270	4.37	0.67
285	4.36	0.62
300	4.42	0.58];

tspan = tgc(:,1);

alpha = 0.09; % min-1
beta = 36; % 10e-9 min-1
h = 5.1615; % mmol/L
K = 710; % 10e-9 

k01 = 0.062;
k21 = 0.053;
k12 = 0.051;
V = 4.18; % L

CP1b = tgc(1,3);

x0(1) = 0;
x0(2) = CP1b;
x0(3) = CP1b*k21/k12;

p = [alpha, beta, h, K, k01, k21, k12];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' Enter to continue with MLE'); disp(' ')
pause
if 1
p_init = [alpha, beta, h, K];
%3 known parameters:
p_fix = [k01,k21, k12, V, x0(1),x0(2),x0(3)];
end

lb = 0*ones(size(p_init));%[0 0 0 0];
ub = [];%ones(size(p_init));%[];
options = optimset('Display','iter','TolFun', 1e-4,...%'iter' default:1e-4
'TolX',1e-5,... %default: 1e-4
'MaxFunEvals' , [],... %800 default:
'LargeScale','on'); %default: on

plt = 0;

[p_est,resnorm,RESIDUAL,exitflag,OUTPUT,LAMBDA,Jacobian]= lsqnonlin(@err_fn,p_init,lb,ub,options,p_fix,tspan,tgc); 
disp(' ')

varp = resnorm*inv(Jacobian'*Jacobian)/length(tspan);
stdp = sqrt(diag(varp)); %The standard deviation is the square root of the variance

% p_est = [alpha, beta, h, K];
% p_fix = [k01,k21, k12, V, x0(1),x0(2),x0(3)];
p = [p_est(1), p_est(2), p_est(3), p_est(4), p_fix(1), p_fix(2),p_fix(3),p_fix(4)];
%x0 = [x0]
x0 = [p_fix(5), p_fix(6), p_fix(7)];

plt =1;
cpep = cpep_sim(tspan,x0,tgc, p,plt);
figure(1); subplot(121)
h1 = plot(tspan,tgc(:,3),'or', 'Linewidth',2);

disp(' Parameters:')
disp([' alpha = ', num2str(p_est(1)), ' +/- ', num2str(stdp(1))])
disp([' beta = ', num2str(p_est(2)), ' +/- ', num2str(stdp(2))])
disp([' h = ', num2str(p_est(3)), ' +/- ', num2str(stdp(3))])
disp([' K = ', num2str(p_est(4)), ' +/- ', num2str(stdp(4))])
disp([' k01 = ', num2str(p_fix(1))])
disp([' k21 = ', num2str(p_fix(2))])
disp([' k12 = ', num2str(p_fix(3))])
disp([' V = ', num2str(p_fix(4))])
disp(' ')
disp(' Initial conditions:')
disp([' Y(0) = ',  num2str(p_fix(5))])
disp([' CP1(0) = ', num2str(p_fix(6))])
disp([' CP2(0) = ', num2str(p_fix(7))])
disp(' ')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function e = err_fn(p_var, p_fix,tspan,tgc)
p = [p_var(1), p_var(2), p_var(3),p_var(4), p_fix(1),p_fix(2),p_fix(3),p_fix(4)];

x0 = [ p_fix(5), p_fix(6), p_fix(7)];
cpep = cpep_sim(tspan,x0,tgc, p,0);
%LSQNONLIN: objective function should return the model error
e = cpep-tgc(:,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cpep = cpep_sim(tspan,x0,tgc, p,plt)

ode_options = [];
[t,x] = ode45(@ode_fn,tspan,x0,ode_options,tgc,p);

cpep = x(:,2);

if plt==1
cpep_plt(tspan,x,tgc,p)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dx = ode_fn(t,x,tgc,p)
%CPEP_ODE ODE's of C-peptide minimal model

% p_est = [alpha, beta, h,K, k01,k21, k12,V];
alpha = p(1);
beta = p(2);
h = p(3);
K = p(4);
k01 = p(5);
k21 = p(6);
k12 = p(7);
V = p(8);
Gb = tgc(1,2);
CP1b = tgc(1,3);

g = interp1(tgc(:,1),tgc(:,2),t);    % using experimental glucose value over testing time

if g >= h
    dx1 = - alpha * ((x(1) - beta*(g-h)));
else
    dx1 =  - alpha * x(1);
end

if diff(g)>0 & g>Gb
    ISRd = K * diff(g)./diff(t);
else
    ISRd = 0;
end

ISRb = k01 * CP1b * V;
ISR = ISRb + ISRd + x(1);

dx2 = ISR -(k01 + k21)* x(2) + k12*x(3);
dx3 = k21 * x(2) - k12 * x(3);

dx = [dx1;dx2;dx3];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cpep_plt(tspan,x,tgc,p)
Gb = tgc(1,2);
CP1b = tgc(1,3);

figure
subplot(121); h = plot(tspan,x(:,2), '-','Linewidth',2); hold on
plot( [tspan(1) tspan(end)], [CP1b CP1b], '--k','Linewidth',1.5)
%baseline level
ylabel('Cpep level [uU/mL]')
legend(h, 'model simulation')

g = interp1(tgc(:,1),tgc(:,2), tspan); %reconstruct used input signal
subplot(122); plot(tspan,g, '--or','Linewidth',2); hold on
plot( [tspan(1) tspan(end)], [Gb Gb], '--k','Linewidth',1.5)
%baseline level
ylabel('glucose level [mg/dL]'); xlabel('time [min]')
legend('interpolated measured test data')
