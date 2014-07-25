function runode1
    tspan = 0:0.1:100;
    init_cond1 = [.1];
    
[t, x] = ode45(@ode_fun1, tspan, init_cond1);

    save('ode1data');

    figure; plot(t, x(:, 1),'b');; 
    
end

function dxdt = ode_fun1(t, x)
    dxdt = zeros(1, 1);
    n=2;
    k=.18; %
    beta=0.36;

    dx1dt = beta * (x(1)^n / (k^n + x(1)^n))- 0.18*x(1);
    dxdt(1) = dx1dt;
end