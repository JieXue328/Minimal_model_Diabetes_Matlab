function [mean,input2] = ModifiedInput_061514(j)
input2 = [];
mean = [];
for i = 1:12
    [X,Y] = MeanWithNan_061514(j(:,i));
    mean = [mean X];
    input2 = [input2 Y];
end