function [m,i] =MeanWithNan_061514(i)
boolean = isnan(i);
NumNotNan = length(i)- sum(boolean);
i(boolean) = 0;
m = sum(i)/ NumNotNan;
i(boolean)= m


