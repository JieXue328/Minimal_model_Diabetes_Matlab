function myFit_070514 % myFit weighed least squares fit

%% create the first half of the data
xdata1 = 0:.01:1;
ydata1 = exp(-3*xdata1) + randn(size(xdata1))*0.05;
weight1 = ones(size(xdata1))*1;

%% create the second half of the data
xdata2 = 1:.01:2;
ydata2 = (xdata2-1).^2 + randn(size(xdata2))*.05;
weight2 = ones(size(xdata2))*10;

%% combine the two data sets
xdata = [xdata1 xdata2];
ydata = [ydata1 ydata2];
weight = [weight1 weight2];

%% call /lsqnonlin/
parameter_hat = lsqnonlin(@mycurve,[.1 -1],[],[],[],xdata,ydata)

%% plot the original data and fitted function
plot(xdata,ydata,'b.')
hold on
fitted = exp(parameter_hat(1).*(parameter_hat(2)+ xdata));
plot(xdata,fitted,'r')
xlabel('x');ylabel('y')
legend('Data','Fit')

%% function that reports the error
function err = mycurve(parameter,real_x,real_y)
fit = exp(parameter(1).*(real_x + parameter(2)));
err = fit - real_y

% weight the error according to the /WEIGHT/ vector
err_weighted = err.* weight;
err = err_weighted
end
end
