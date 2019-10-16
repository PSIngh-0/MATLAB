% Parminder Singh
% 004834219

function [x, y] = rand_walk(x, y, bound)
% This function has the inputs x, y, and boud with the output of [x, y]
% rand_walk is a function which inputs position information and randomly
% decides whether a particle should move up, down, left, right, or remain
% still.
% The porbability for each of these is equal
% If the particle is at a boundary its position will remain unchanged

% Randomly pick a number between 0 and 1 but not including 0 or 1
n = rand;
p = 0.2;
pStill = 1 - 4*p;


if n < pStill % Particle remains still
    x = x + 0;
    y = y + 0;
    
    if x > bound(1) % Boundary conditions
        x = bound(1);
    end
    
    if x < bound(2)
        x = bound(2);
    end
    
    if y > bound(3)
        y = bound(3);
    end
    
    if y < bound(4)
        y = bound(4);
    end 
    
elseif pStill < n && n < 2*p % Particle moves up
    y = y + 1;
    
    if y > bound(3) % Boundary conditions
        y = bound(3);
    end
    
elseif 2*p < n && n < 3*p % Particle moves right
    x = x + 1;
    
    if x > bound(1) % Boundary conditions
        x = bound(1);
    end
    
elseif 3*p < n && n < 4*p % Particle moves down
    y = y - 1;
    
    if y < bound(4) % Boundary conditions
        y = bound(4);
    end
    
elseif 4*p < n && n < 5*p  % Particle moves left
    x = x - 1;
    
    if x < bound(2) % Boundary conditions
        x = bound(2);
    end
    
end

end