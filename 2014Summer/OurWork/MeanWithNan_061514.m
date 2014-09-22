function [m,i] =MeanWithNan_061514(i)   
boolean = isnan(i);                        % find the NaN values
NumNotNan = length(i)- sum(boolean);       % number of non-NaN values per row/column
i(boolean) = 0;                            % replace NaN with 0
m = sum(i)/ NumNotNan;                     % calculate & return the mean value per row/column
i(boolean)= m                              % replace initial NaN to mean value $ return the modified vector


