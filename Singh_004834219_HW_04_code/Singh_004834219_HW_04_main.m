%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Homework 04: Main Matlab Script
% [1] The Split-and-Average Problem
% [2] Runge-Kutta Radioactivity
%
% Author: Parminder Singh
% Date:   07/20/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clear Cache
clear all
close all
clc

%% Homework 04 Problems List
fprintf('Homework 04 Problems:\n\n')
fprintf('1: The Split-and-Average Problem\n')
fprintf('2: Runge-Kutta Radioactivity\n\n')

problem = input('Please enter the problem number: \n');
    while (problem <= 0 || problem > 2) %error check for problem input
        fprintf('\nERROR: Invalid problem number! Please try again.\n\n');
        fprintf('1: The Split-and-Average Problem\n')
        fprintf('2: Runge-Kutta Radioactivity\n\n')
        problem = input('Please enter the problem number: \n');
    end
switch(problem)

    
%% Problem Number 1: The Split-and-Average Problem ========================
case 1

% Initial array
x = [1, 0, 1, 0];
y = [0, 1, 1, 0];
w = [2, 2, 2];

x1 = x; % Used later for original value plot
y1 = y;

Max = 2; % In order to initialize the while loop

while Max >= (1*10^(-3)) 
    xs = splitPts(x); 
    ys = splitPts(y);
    
    xa = averagePts(xs, w);
    ya = averagePts(ys, w);
    
    x = xa;
    y = ya;
    
    dx = xa - xs;
    dy = ya - ys;
    
    Max = max(sqrt(dx.^2 + dy.^2)); 
end

%Plots
figure(1)
        hold on
        plot(x1, y1, 'co', x, y, 'ro')
        legend('Initial Distribution', 'Final Distribution')
        title('The Split Average Problem')

%% Problem Number 2: Runge-Kutta Radioactivity ============================
case 2

% Time values
dt = [1, 0.1, 0.01];
tFinal = 15;
y = 1; % set in order to initiate for loop

% Begin table
fprintf('   dt\t')
fprintf('      RK1\t')
fprintf('     RK2\t')
fprintf('     RK4\t\n')

% for loop to initialize plot values for each method
for k = 1:1:3
    timeStep = dt(k);
    t = 0:timeStep:tFinal;
    steps = length(t) - 1;
        
    exact = zeros(1, steps + 1);
    RK1 = zeros(1, steps + 1);    
    RK2 = zeros(1, steps + 1);    
    RK4 = zeros(1, steps + 1);    
    
    exact(1) = y;
    RK1(1) = y;                   
    RK2(1) = y;
    RK4(1) = y;
    
% for loop for each method
    for j = 1:1:steps
        RK1(j+1) = advanceRK(RK1(j), timeStep, 1);
        RK2(j+1) = advanceRK(RK2(j), timeStep, 2);
        RK4(j+1) = advanceRK(RK4(j), timeStep, 4);
        
       
        timeExact = j * timeStep;
        exact(j+1) = y * exp(f(timeExact));
    end
    
% Error and 
    errRK1 = RK1 - exact;
    errRK2 = RK2 - exact;
    errRK4 = RK4 - exact;
    
    avgRK1 = mean(errRK1);
    avgRK2 = mean(errRK2);
    avgRK4 = mean(errRK4);
    
    fprintf('%.2f:\t', timeStep)
    fprintf('%.2i\t', avgRK1)
    fprintf('%.2i\t', avgRK2)
    fprintf('%.2i\t\n', avgRK4)
    
% Plots for all three methods
    figure                     
    plot (t, exact, t, RK1, t, RK2, t, RK4, 'linewidth', 4, 'linewidth', 2, 'linewidth', 2, 'linewidth', 2)
    legend ('Exact','Runge-Kutta 1','Runge-Kutta 2','Runge-Kutta 4')
    
% Title of  plots given using the first for loop
    switch k
        case 1
            title('dt = 1 s')
        case 2
            title('dt = 0.1 s')
        case 3
            title('dt =  0.01 s')
    end
    
end

%%=========================================================================
end
