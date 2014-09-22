function Ins_mm_mle_A_072514

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Input
% glucose concentration in plasma [mg/dL];
%time (minutes) glucose level (mg/dl)
tgi = [2 233.40	55.96
3	238.30	177.00
4	239.19	325.00
5	224.39	340.01
6	222.14	311.56
8	214.41	266.56
10	205.73	224.14
12	200.58	211.73
14	193.50	218.47
16	187.71	210.21
19	179.99	205.05
22	173.55	193.17
25	165.19	181.29
30	156.18	184.43
40	137.84	169.48
50	120.95	157.13
60	104.54	139.60
70	100.68	132.42
80	93.92	123.17
90	85.71	106.15
100	83.78	102.59
110	80.89	93.86
120	80.41	94.96
140	78.96	94.05
160	79.92	86.93
180	81.37	83.96
210	83.30	75.35
240	83.78	74.51];

%Fixed Ins_Mini_Model initial conditions & model parameters
h =93;
n = 0.18;
gamma1 = 0.0076;
Gb = tgi(1,2);
Ib = tgi(1,3);
I0 = 50;

p = [gamma1, h, n,I0];

x0 = I0;

tspan = tgi(:,1);
tgi = [tgi(:,[1:2]) (tgi(:,3)-tgi(1,3))];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' Enter to continue with MLE'); disp(' ')
pause
if 1
p_init = [gamma1, h, n,I0];
%3 known parameters:
p_fix = [];
end


lb = 0*ones(size(p_init));%[0 0 0 0];
ub = [];%ones(size(p_init));%[];
options = optimset('Display','iter','TolFun', 1e-4,...%'iter' default:1e-4
'TolX',1e-5,... %default: 1e-4
'MaxFunEvals' , [],... %800 default:
'LargeScale','on'); %default: on

plt = 0;

[p_est,resnorm,RESIDUAL,exitflag,OUTPUT,LAMBDA,Jacobian]= lsqnonlin(@err_fn,p_init,lb,ub,options,p_fix,tspan,tgi); 
disp(' ')
%end
%Accuracy:
%lsqnonlin returns the jacobian as a sparse matrix
varp = resnorm*inv(Jacobian'*Jacobian)/length(tspan);
stdp = sqrt(diag(varp)); %The standard deviation is the square root of the variance

% p_est = [gamma1, h, n,I0]
p = [p_est(1), p_est(2), p_est(3), p_est(4)];
%x0 = [x0]
x0 = [p_est(4)];

plt =1;

ins = ins_sim(tspan,x0,tgi, p,plt);
figure(1); subplot(121)
h1 = plot(tspan,tgi(:,3),'or', 'Linewidth',2);

disp(' Parameters:')
disp([' gamma1 = ', num2str(p_est(1)), ' +/- ', num2str(stdp(1))])
disp([' h = ', num2str(p_est(2)), ' +/- ', num2str(stdp(2))])
disp([' n = ', num2str(p_est(3)), ' +/- ', num2str(stdp(3))])
disp([' I0 = ', num2str(p_est(4)), ' +/- ', num2str(stdp(4))])

disp(' Initial conditions:')
disp([' x0 = ', num2str(p_est(4)), ' +/- ', num2str(stdp(4))])
disp(' ')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function e = err_fn(p_var, p_fix,tspan,tgi)
p = [p_var(1), p_var(2), p_var(3),p_var(4)];

x0 = [p_var(4)];
ins = ins_sim(tspan,x0,tgi, p,0);

weight1 = (tspan<5)*0;
weight2 = (tspan >=5 & tspan<=120)*10;
weight3 = (tspan>120)*1;
weight = weight1 + weight2 + weight3;

e = ins-tgi(:,3);

err_weighted = e .* weight;
e = err_weighted;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ins = ins_sim(tspan,x0,tgi, p,plt)

ode_options = [];
[t,x] = ode45(@ode_fn,tspan,x0,ode_options,tgi,p);

ins = x;

if plt==1
ins_plt(tspan,x,tgi,p)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dx = ode_fn(t,x,tgi,p)
%INS_ODE ODE's of insulin minimal model

gamma1 = p(1);
h = p(2);
n = p(3);
I0 = p(4);

g = interp1(tgi(:,1),tgi(:,2),t);  % using experimental glucose value over testing time

if g > h
    dx = gamma1 * (g - h)*t - (n * x);
else
    dx =  - (n * x);
end
           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ins_plt(tspan,x,tgi,p)
Gb = tgi(1,2);
Ib = tgi(1,3);

figure
subplot(121); h = plot(tspan,x, '-','Linewidth',2); hold on
plot( [tspan(1) tspan(end)], [Ib Ib], '--k','Linewidth',1.5)
%baseline level
ylabel('insulin level [uU/mL]')
legend(h, 'Insulin minimal model simulation')

g = interp1(tgi(:,1),tgi(:,2), tspan); %reconstruct used input signal
subplot(122); plot(tspan,g, '--or','Linewidth',2); hold on
plot( [tspan(1) tspan(end)], [Gb Gb], '--k','Linewidth',1.5)
%baseline level
ylabel('glucose level [mg/dL]'); xlabel('time [min]')
legend('interpolated measured test data')
%figure




