% Parminder Singh
% 004834219

function [y] = advanceRK(y, dt, method)
% advanceRK uses three different Runge-Kutta methods to advance a 
% discretized solution by one time stamp

switch (method)
    case 1 
        c1 = dt * f(y);
        y = y + c1;
    case 2 
        c1 = dt * f(y);
        c2 = dt * f(y + (0.5 * c1));
        y = y + c2;
    case 4 
        c1 = dt * f(y);
        c2 = dt * f(y + (0.5 * c1));
        c3 = dt * f(y + (0.5 * c2));
        c4 = dt * f(y + c3);
        y = y + ((1/6) * c1) + ((1/3) * c2) + ((1/3) * c3) + ((1/6) * c4);
end

end

