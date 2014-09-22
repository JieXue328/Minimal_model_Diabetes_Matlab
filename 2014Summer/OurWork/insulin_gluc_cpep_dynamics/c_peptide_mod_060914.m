function c_peptide_mod_060914

%Fixed initial conditions & model parameters

x0(1) = 0.568;
x0(2) = 1.11;
Ib = 0.057;
Cb = 0.379;
Rb = 8.9 
f = 0.818; 
kI = 0.122; 
kc = 0.028; 
r = 1.28

p = [f, kI, kc];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% call ODE function
tspan = 0:5:180;
ode_options = [];
[t,x] = ode45(@cpeptide_ode,tspan,x0,ode_options,p);

insul = x(:,1);  
c_pept = x (:,2);    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plotting
figure
fig1 = plot(tspan,insul,'r', 'Linewidth',2); hold on   
plot( [tspan(1) tspan(end)], [Ib Ib], '--k','Linewidth',1.5)    
xlabel('time [min]'); ylabel('insulin level [pmol/mL]')
title('INSULIN_TIME_PROFILE')

figure
fig2 = plot(tspan,c_pept,'g', 'Linewidth',2); hold on   
plot( [tspan(1) tspan(end)], [Cb Cb], '--k','Linewidth',1.5)    

xlabel('time [min]'); ylabel('C-peptide level [pmol/mL]')
title('C-PEPTIDE_TIME_PROFILE')

figure 
fig3 = plot(tspan, r,'b', 'Linewidth',2)                     
xlabel('time [min]'); ylabel('insulin secrete rate [1/min]')
title('INSULIN SCRETE RATE_TIME_PROFILE')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dxout = cpeptide_ode(t,xin,p)

f = p(1);
kI = p(2);
kc = p(3);

%ode's
dI = f* 1.28 - kI*xin(1);          % xin(1) denotes I(t); xin(2) denotes C(t); xin(3) denotes r(t) 
dC = 1.28  - kc* xin(2); 
dxout = [dI; dC];
