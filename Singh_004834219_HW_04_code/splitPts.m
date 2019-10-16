% Parminder Singh
% 004834219

function [xs] = splitPts(x)
% The function splitPts splits each point in an array by creating an empty
% array and filling a new array with an extra spot between each number that
% was in the old array.
% This doubles the size of the original array

% Arrays
xa = length(x);
n = 2 * xa;
xs = zeros(1,n);

% Split average calulcation within array
for k = 1:1:xa
    newXA = (k * 2) - 1; % making new array with a spot between each number
    xs(newXA) = x(k); 
    if k == xa % if the number is the last number
        split = (x(k) + x(1))/2; 
    else 
        split = (x(k) + x(k + 1))/2; 
    end
    xs(newXA + 1) = split; % creating new array with split numbers
end
end

