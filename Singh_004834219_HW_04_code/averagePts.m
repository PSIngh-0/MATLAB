% Parminder Singh
% 004834219

function [xa] = averagePts(xs, w)
% The function averagePts takes a weighted average of all adjacent pairs
% in an array.
% A weight vector should be defined in the main script so it could be 
% used here. 



if sum(w) == 0
    error('The sum of the weighted averages cannot equate to 0!')
end

w = w/sum(w);
n = length(xs);
xa = zeros(1, n);

for k = 1:1:n       
    if k == 1 %First value
        xa(k) = (w(1) * xs(n)) + (w(2) * xs(k)) + (w(3) * xs(k+1));
    elseif k == n % Last value
        xa(k) = (w(1) * xs(k-1)) + (w(2) * xs(k)) + (w(3) * xs(1));
    else % Any other value 
        xa(k) = (w(1) * xs(k-1)) + (w(2) * xs(k)) + (w(3) * xs(k+1));
    end
end
end

