function Cpep_Deconvol_Model_062414

K1 = 0.057;
K2 = 0.054;
K3 = 0.06;
Dose = 4;
Vp = 65.5;

alpha = (K1 + K2 + K3 + sqrt((K1 + K2 + K3)^2 - 4*K2*K3))/2;
beta = (K1 + K2 + K3 - sqrt((K1 + K2 + K3)^2 - 4*K2*K3))/2;

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
cpep = (Dose * ((alpha-K2)*exp(-alpha*tspan)-(beta-K2)*exp(-beta*tspan)))/ (Vp * (alpha - beta))

figure                % plot cpep-time profile
h = plot(tspan,cpep,'-', 'Linewidth',2);
legend(h,'simulated cpeptide-time profile');
xlabel('time'); ylabel('Cpeptide level');
title('CPEPTIDE DECONVOLUTION MODEL')