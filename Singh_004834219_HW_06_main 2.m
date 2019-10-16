%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Homework 06: Main Matlab Script
% [1] The Game of Life
% [2] Euler-Bernoulli Beam Bending
%
% Author: Parminder Singh
% Date:   08/05/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clear Cache
clear all
close all
clc

%% Homework 04 Problems List
fprintf('Homework 06 Problems:\n\n')
fprintf('1: The Game of Life\n')
fprintf('2: Euler-Bernoulli Beam Bending\n\n')

problem = input('Please enter the problem number: \n');
    while (problem <= 0 || problem > 2) %error check for problem input
        fprintf('\nERROR: Invalid problem number! Please try again.\n\n');
        fprintf('1: The Game of Life\n')
        fprintf('2: Euler-Bernoulli Beam Bending\n\n')
        problem = input('Please enter the problem number: \n');
    end
switch(problem)

    
%% Problem Number 1: The Game of Life =====================================
case 1
    
rng('shuffle');    

x_grid = 200;
y_grid = 150;
n = 300;
time = zeros(1, n);

grid = rand(y_grid, x_grid);

for j = 1:1:y_grid
    for i = 1:1:x_grid
        
        if grid(j, i) < 0.1
            grid(j, i) = 1;
            
        else
            grid(j, i) = 0;
        end
    end
end

for k = 1:1:n
    new_grid = zeros(j, i);
    for j = 1:1:y_grid
        for i = 1:1:x_grid
            N = j - 1;
            E = i + 1;
            S = j + 1;
            W = i - 1;
            if N < 1
                N = 150;
            end
            if S > y_grid
                S = 1;
            end
            if W < 1
                W = x_grid;
            end
            if E > x_grid
                E = 1;
            end
            neighbors = grid(N, E) + grid(N, W) + grid(S, E)...
                      + grid(S, W) + grid(N, i) + grid(S, i)... 
                      + grid(j, E) + grid(j, W);
            if grid(j, i) == 1 && (neighbors == 2 || neighbors == 3)
                new_grid(j, i) = 1;
            end 
            if grid(j, i) == 1 && (neighbors < 2 || neighbors > 3)
                new_grid(j, i) = 0;
            end 
            if grid(j, i) == 0 && neighbors == 3
                new_grid(j,i) = 1; 
            end
        end
    end

% Data for time vs population plot
sum_grid(k) = sum(sum(grid));
grid = new_grid;
time_steps = linspace(1, 300, 300);


% Animation Plot
figure(1)
imagesc(grid)
title('The Game of Life over 300 Generations')
set(gcf, 'Position', [30 350 850 450])
set(gca, 'LineWidth', 3, 'FontSize', 20)
axis equal
drawnow

end 

% Plot for time vs population
plot(time_steps, sum_grid);
xlabel('Time');
ylabel('Population');
title('Population vs. Time');

%% Problem Number 2: Euler-Bernoulli Beam Bending =========================
case 2

% Variables
nodes = 20;
L = 1;
P = -2000;
r2 = 0.013;
r1 = 0.011;
d = 0.75;

% Equations
E = 7e10;
I = (pi/4)*(r2^4 - r1^4);

% Create matrix A
A = zeros(nodes, nodes);

% Boundary conditions
A(1, 1) = 1;
A(nodes, nodes) = 1;

% For loop for assigning matrix values to A
for k = 2:1:(nodes - 1)
    A(k, k-1) = 1;
    A(k, k) = -2;
    A(k, k+1) = 1;
end

% Vectors 
x = linspace(0, 1, nodes);
dx = x(2) - x(1);

% Create matrix B
B = zeros(nodes, 1);

% For loop for assigning matrix values to B
for j = 1:1:nodes
    if x(j) <= d
        M = ((-P)*(L - d)*x(j))/L;
    else
        M = ((-P)*d*(L - x(j)))/L;
    end
    B(j) = ((dx^2)*M)/(E*I);
end

% Beam bending error calculation
y = A\B;
y_max_experimental = min(y);
c = min(d,(L - d));
y_max_theoretical = (P*c*((L^2 - c^2)^1.5))/(9*(sqrt(3))*E*I*L);
error = abs(y_max_experimental - y_max_theoretical)*100;
fprintf('\nMaximum displacement error: %.10f \n\n',error)

% Plot for beam bending
plot(x, y,'o-')
xlabel('Distance along Beam (m)')
ylabel('Deflection of Beam (m)')
title('Deflection of a Beam of Length L Under a Single Point Load at d = 0.75L')

%%=========================================================================
end
