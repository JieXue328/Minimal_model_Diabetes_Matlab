function Ins_Cpep_Mini_Model_B_062414

%Fixed Ins_Mini_Model initial conditions & model parameters
h = 4.6;
m1 = 0.64;
alpha1 = 0.09;
beta1 = 57;
y0(1) = 701;
y0(2) = 0;
y0(3) = 0;
n = 0.19;
V1 = 9.2
p = [h,m1,alpha1,beta1,n,V1];

%Fixed Cpep_Mini_Model initial conditions & model parameters
h = 4.6
m2 = 0.68;
alpha2 = 0.14;
beta2 = 8.2;
z0(1) = 1863;
z0(2) = 0;
z0(3) = 0;
z0(4) = 0;
k01 = 0.062;
k21 = 0.053;
k12 = 0.051;
q = [h, m2, alpha2, beta2, k01,k21,k12]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Input
% glucose concentration in plasma [mg/dL];
%time (minutes) glucose level (mg/dl)
tg = [ 0 92 
    2 350 
    4 287 
    6 251 
    8 240 
    10 216
    12 211 
    14 205
    16 196 
    19 192 
    22 172 
    27 163 
    32 142 
    42 124 
    52 105 
    62 92 
    72 84 
    82 77 
    92 82 
    102 81 
    122 82 
    142 82 
    162 85 
    182 90 ];
tspan = tg(:,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ode_options = [];
[t,y] = ode45(@ins_ode,tspan,y0,ode_options,tg,p);
[t,z] = ode45(@cpep_ode,tspan,z0,ode_options,tg,q);

ins = y(:,3);
cpep = z(:,3);

figure         % plot insulin-time profile
h1 = plot(tspan,ins,'-g', 'Linewidth',2);
legend(h1,'simulated insulin-time profile');
xlabel('time'); ylabel('insulin level');
title('INSULIN MINMAL MODEL')

figure          % plot cpep-time profile
h2 = plot(tspan,cpep,'-g', 'Linewidth',2);
legend(h2,'simulated cpeptide-time profile');
xlabel('time'); ylabel('Cpeptide level');
title('CPEPTIDE MINMAL MODEL')

figure          % plot insulin secretion rate-time profile
isr = m2 * z(:,1);
h3 = plot(tspan,isr,'-g', 'Linewidth',2);
xlabel('time'); ylabel('insulin secretion rate');
title('INSULIN SECRETION RATE TIME PROFILE')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dy = ins_ode(t,y,tg,p)
%INS_ODE ODE's of insulin minimal model

h = p(1);
m1 = p(2);
alpha1 = p(3);
beta1 = p(4);
n = p(5);
V1 = p(6);

g = interp1(tg(:,1),tg(:,2),t);   % using experimental glucose value over testing time

dX = -m1*y(1) + y(2)
dY = -alpha1 * (y(2)-beta1*(g-h))
dI = -n*y(3) + m1*y(1)          % plasma insulin secretion and kinetics model
dy = [dX;dY;dI]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dz = cpep_ode(t,z,tg,q)
%CPEP_ODE ODE's of C-peptide minimal model

h = q(1);
m2 = q(2);
alpha2 = q(3);
beta2 = q(4);
k01 = q(5);
k21 = q(6);
k12 = q(7);

g = interp1(tg(:,1),tg(:,2),t);    % using experimental glucose value over testing time

dX = -m2*z(1) + z(2)
dY = -alpha2 * (z(2)-beta2*(g-h))
dCP1 = -(k01 + k21)* z(3) + k12*z(4) + m2*z(1); % plasma c-peptide secretion and kinetics model
dCP2 = k21 * z(3) - k12*z(4);
dz = [dX;dY;dCP1;dCP2];
