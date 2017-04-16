function [mean,input2] = ModifiedInput_061514(j)
input2 = [];                      % initiaize modified matrix
mean = [];                        % initialize the mean vector
for i = 1:length(j(1,:));                      % ovall 12 columns
    [X,Y] = MeanWithNan_061514(j(:,i)); % calucate & return mean value and modified vector for each column
    mean = [mean X];                    % generate a row vector of mean values
    input2 = [input2 Y];                % generate modified matrix
end