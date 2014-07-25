function compartment_072114
data =[0	92.78	25.09
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

t = data(:,1);
glu = data(:,2);
ins = data(:,3);

figure; plot(t,glu,'.b',t,ins,'.r');hold on
xlabel('Time [s]'); ylabel(' [mM]')

N = length(ins);

%% simulate model with initial parameter values
Ktrans = 2e-3; ve = 0.1;
p = [Ktrans,ve];
[y,t] = compart_lsim(t,glu,p);
plot(t,y,'k')

%% estimate parameters and variances
optim_options = optimset('Display','iter','TolFun',1e-6,'TolX',1e-6,'LevenbergMarquardt','on');
p0 = [Ktrans,ve];
[p,resnorm,residual,exitflag,OUTPUT,LAMBDA,Jacobian] = lsqnonlin(@compart_err,p0,[],[],optim_options,t,glu,ins);
disp(' ')
p

% estimate parameter variance
Jacobian = full(Jacobian);
varp = resnorm*inv(Jacobian'*Jacobian)/N;
stdp = sqrt(diag(varp));
stdp = 100*stdp'./p;
disp([' Ktrans: ', num2str(p(1)),'+/-',num2str(stdp(1)),'%'])
disp([' ve: ', num2str(p(2)),'+/-',num2str(stdp(2)),'%'])


%% simulate estimated model
[y,t] = compart_lsim(t,glu,p);
figure;subplot(211);plot(t,ins,'.r',t,y,'b')
xlabel('Time [s]'); ylabel(' [mM]');
legend('data','model')
xi = ins(:)-y(:);
subplot(212);plot(t,xi)
xlabel('Time [s]');legend('residuals \xi')
assen=axis;

%% function to simulate compartment_model
function [y,t] = compart_lsim(t,glu,p);
Ktrans = p(1); ve = p(2);
num = [Ktrans*ve];
den = [ve,Ktrans];
sys = tf(num,den);
[y,t] = lsim(sys,glu,t);
end

%% function to calculate MLE error
    function xi =compart_err(p,t,glu,ins)
        [y,t]= compart_lsim(t,glu,p);
        xi = ins(:)-y(:);
    end
end







