function gluc_oral_mm_mle_081414

%time (minutes) glucose level (mmol/L) insulin level (pmol/L)
tgi = [0	5.39	94.72
10	5.93	180.27
20	6.95	362.36
30	8.05	480.69
45	9.04	586.25
60	9.44	617.12
75	9.27	735.43
90	8.84	602.37
105	8.25	722.5
120	7.58	673.23
135	7.05	541.99
150	7.03	642.08
165	6.54	447.09
180	6.01	348.64
195	5.45	257.47
210	5.04	189.99
225	4.78	182.61
240	4.48	104.19
255	4.34	104.1
270	4.37	113.12
285	4.36	82.06
300	4.42	74.68];

glu_exp = tgi(:,2);
ins_exp = tgi(:,3);
tspan = tgi(:,1);
Gb = tgi(1,2);
Ib = tgi(1,3);

x0(1) = Gb;
x0(2) = 0;
p1 = 2.6e-2;
p2 = 0.025;
p3 = 5.0e-4;
V = 1.5;

alpha = rand(1,7);
alpha= [0 alpha];
p = [p1 p2 p3 V alpha];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' Enter to continue with MLE'); disp(' ')
pause
if 1
%11 unknown model parameters:
p_init = [p1 p2 p3 V alpha(2:end)];
%5 known parameters:
p_fix = [Gb Ib alpha(1) x0(1) x0(2)];
end

lb = 0*ones(size(p_init));%[0 0 0 0];
ub = [];%ones(size(p_init));%[];
options = optimset('Display','iter','TolFun', 1e-4,...%'iter' default:1e-4
'TolX',1e-5,... %default: 1e-4
'MaxFunEvals' , [],... %800 default:
'LargeScale','on'); %default: on
plt = 0;

[p_est,resnorm,RESIDUAL,exitflag,OUTPUT,LAMBDA,Jacobian]= lsqnonlin(@err_fn,p_init,lb,ub,options,p_fix,tspan,tgi, plt); 
disp(' ')

varp = resnorm*inv(Jacobian'*Jacobian)/length(tspan);
stdp = sqrt(diag(varp)); %The standard deviation is the square root of the variance
%p_init = [p1 p2 p3 V alpha(2:end)];
%p_fix = [Gb Ib alpha(1) x0(1) x0(2)];
p = [p_est(1), p_est(2), p_est(3), p_est(4),p_fix(1),p_fix(2), p_fix(3),p_est(5:end)];
%x0 = [G0, X0]
x0 = [p_fix(4), p_fix(5)];
disp(' Parameters:')
disp([' p1 = ', num2str(p_est(1)), ' +/- ', num2str(stdp(1))])
disp([' p2 = ', num2str(p_est(2)), ' +/- ', num2str(stdp(2))])
disp([' p3 = ', num2str(p_est(3)), ' +/- ', num2str(stdp(3))])
disp([' V = ', num2str(p_est(4)), ' +/- ', num2str(stdp(4))])
disp([' Gb = ', num2str(p_fix(1))])
disp([' Ib = ', num2str(p_fix(2))])
disp([' alpha0 = ', num2str(p_fix(3))])
disp([' alpha1 = ', num2str(p_est(5)), ' +/- ', num2str(stdp(5))])
disp([' alpha2 = ', num2str(p_est(6)), ' +/- ', num2str(stdp(6))])
disp([' alpha3 = ', num2str(p_est(7)), ' +/- ', num2str(stdp(7))])
disp([' alpha4 = ', num2str(p_est(8)), ' +/- ', num2str(stdp(8))])
disp([' alpha5 = ', num2str(p_est(9)), ' +/- ', num2str(stdp(9))])
disp([' alpha6 = ', num2str(p_est(10)), ' +/- ', num2str(stdp(10))])
disp([' alpha7 = ', num2str(p_est(11)), ' +/- ', num2str(stdp(11))])
disp(' ')
disp(' Initial conditions:')
disp([' G0 = ', num2str(p_fix(4))])
disp([' X0 = ', num2str(p_fix(5))]); disp(' ')
plt = 1;
glu = glu_sim(tspan,x0,tgi, p,plt);
%compare model output with measured data
figure(1); subplot(221)
h1 = plot(tspan,glu_exp,'or', 'Linewidth',2);
%legend(h1, 'experimental test data')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function e = err_fn(p_var, p_fix,tspan,tgi,plt)
%p_init = [p1 p2 p3 V alpha(2:end)];
%p_fix = [Gb Ib alpha(1) x0(1) x0(2)];

p = [p_var(1), p_var(2), p_var(3), p_var(4),p_fix(1),p_fix(2), p_fix(3),p_var(5:end)];
%x0 = [G0, X0]
x0 = [p_fix(4), p_fix(5)];
glu = glu_sim(tspan,x0,tgi, p, 0);
%LSQNONLIN: objective function should return the model error
e = glu-tgi(:,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function glu = glu_sim(tspan,x0,tgi, p,plt)

ode_options = [];
[t,x] = ode45(@glu_ode,tspan,x0,ode_options, tgi, p);

%Output
glu = x(:,1);
if plt==1
glu_plt(tspan,x,tgi,p)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dxout = glu_ode(t,xin,tgi, p)

%p = [p_est(1), p_est(2), p_est(3), p_est(4),p_fix(1),p_fix(2), p_fix(3),p_est(5:end)];

p1 = p(1); 
p2 = p(2); 
p3 = p(3); 
V = p(4); 
Gb = p(5); 
Ib = p(6);
alpha = p(7:end);

u = interp1(tgi(:,1),tgi(:,3), t);

t_ra = [0,10,30,60,90,120,180,300]';

%ode's
for i = 1:7
    if t >= t_ra(i) & t <= t_ra(i+1)
        Ra = alpha(i) + (alpha(i+1)-alpha(i))*(t-t_ra(i))/(t_ra(i+1)-t_ra(i));
    end
end
dG = -(p1 + xin(2))* xin(1) + p1 * Gb + Ra/V;
dX = -p2*xin(2) + p3*(u-Ib);
dxout = [dG; dX];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function glu_plt(tspan,x,tgi,p)
Gb = p(5);
Ib = p(6);

figure
subplot(221); h = plot(tspan,x(:,1), '-','Linewidth',2); hold on
plot( [tspan(1) tspan(end)], [Gb Gb], '--k','Linewidth',1.5)
%baseline level
ylabel('glucose level [mmol/L]')
legend(h, 'model simulation')

u = interp1(tgi(:,1),tgi(:,3), tspan); %reconstruct used input signal
subplot(222); plot(tspan,u, '--or','Linewidth',2); hold on
plot( [tspan(1) tspan(end)], [Ib Ib], '--k','Linewidth',1.5)
%baseline level
ylabel('insulin level [pmol/L]'); xlabel('time [min]')
legend('interpolated measured test data')
%figure
subplot(223); plot(tspan, x(:,2), 'Linewidth',2)
xlabel('time [min]'); ylabel('interstitial insulin [1/min]')