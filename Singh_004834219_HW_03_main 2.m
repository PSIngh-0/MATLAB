%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Homework 03: Main Matlab Script
%   [1] The Pendulum Physics Problem
%   [2] DNA Analysis
%
%   Author: Parminder Singh
%   Date:   07/15/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clear Cache
clear all
close all
clc

%% Homework 01 Problems List
fprintf('Homework 03 Problems:\n\n')
fprintf('1: The Pendulum Physics Problem\n')
fprintf('2: DNA Analysis\n\n')

problem = input('Please enter the problem number: \n');
    while (problem <= 0 || problem > 2) %error check for problem input
        fprintf('\nERROR: Invalid problem number! Please try again.\n\n');
        fprintf('1: Neighbor Identification\n')
        fprintf('2: The Three Species Problem\n\n')
        problem = input('Please enter the problem number: \n');
    end
switch(problem)

    
%% Problem Number 1: Neighbor Identification ==============================
case 1
% Variables
L = 1;
g = 9.81;
dt = 0.005;
t_f = 20;
t_step = ceil(t_f/dt) + 1;

% Arrays
t = linspace(0, t_f , t_step);

% Explicit Array
theta_0 = zeros(1, t_step);
w_0 = zeros(1, t_step);
a_0 = zeros(1, t_step);
h_0 = zeros(1, t_step);
energy_0 = zeros(1, t_step);
theta_0(1) = pi/3;
w_0(1) = 0;

% Implicit Array
w_1 = zeros(1, t_step);
theta_1 = zeros(1, t_step);
a_1 = zeros(1, t_step);
h_1 = zeros(1, t_step);
energy_1 = zeros(1, t_step);
theta_1(1) = pi/3;
w_1(1) = 0;

% Assign values to arrays
for j = 1 : t_step - 1  
    % Explicit
    w_0(j + 1) = -dt * (g/L) * sin(theta_0(j)) + w_0(j);
    theta_0(j + 1) = dt * w_0(j) + theta_0(j);
    a_0(j) = -(g/L) * sin(theta_0(j));
    h_0(j) = L - cos(theta_0(j)) * L ;
    energy_0(j) = g * h_0(j) + (.5 * (L * w_0(j))^2);
    w_0(j) = w_0(j + 1);
    theta_0(j) = theta_0(j + 1);
    
    % Implicit
    theta_1(j + 1) = -dt^2 * (g/L) * sin(theta_1(j)) + dt * w_1(j) + theta_1(j);
    w_1(j + 1) = (theta_1(j+1) - theta_1(j))/dt;
    a_1(j) = -(g/L) * sin(theta_1(j));
    h_1(j) = L - cos(theta_1(j)) * L ;
    energy_1(j) = g * h_1(j) + (.5 * (L * w_1(j))^2);
    w_1(j) = w_1(j+1);
    theta_1(j) = theta_1(j+1);
end

% Explicit Plot
figure(1)
plot(t, theta_0 , t, w_0, t, a_0);
xlabel('Time (s)');
ylabel('Position (rad), Velocity (rad/s) , Acceleration (rad/s^2)');
legend('Angular Position (rad)','Angular Velocity (rad/s)','Angular Acceleration (rad/s^2)')
title('Explicit Angular Position, Velocity, and Acceleration');

% Explicit Energy vs Time
figure(2)
plot(t ,energy_0);
xlabel('Time (s)');
ylabel('Energy (J)')
title('Explicit Energy over Time');

% Implicit Plot
figure(3)
plot(t, theta_1 , t, w_1, t, a_1);
xlabel('Time (s)');
ylabel('Position (rad), Velocity (rad/s) , Acceleration (rad/s^2)');
legend('Angular Position (rad)','Angular Velocity (rad/s)','Angular Acceleration (rad/s^2)')
title('Implicit Angular Position, Velocity, and Acceleration');

% Implicit Energy vs Time
figure(4)
plot(t , energy_1);
xlabel('Time (s)');
ylabel('Energy (J)')
title('Implicit Energy over Time');

%% Problem Number 2: The Three Speicies Problem ===========================
case 2
% Load DNA array
load('chr1_sect.mat');

% DNA vector length caaulculation and starting point
numBase = length(dna);
start = 0;

% Protein length 
length = zeros(1, numBase);
TAA = 0;
TAG = 0;
TGA = 0;

% Start and stop codons 
for j = 1 : 3 : numBase - 2
    if start == 0
        if dna(j) == 1 && dna (j + 1) == 4 && dna(j + 2) == 3
            start = j;
        end
    elseif (dna(j) == 4) && ((dna(j + 1) == 1 && dna(j + 2) == 1) || (dna(j + 1) == 1 && dna(j + 2) == 3) || (dna(j + 1) == 3 && dna( j + 2) == 1))
            stop = j;
            if (dna(j) == 4) && (dna(j + 1) == 1 && dna(j + 2) == 1)
                TAA = TAA + 1;
            end
            if (dna(j) == 4) && (dna(j + 1) == 1 && dna(j + 2) == 3)
                TAG = TAG + 1;
            end
            if (dna(j) == 4) && (dna(j + 1) == 3 && dna(j + 2) == 1)
                TGA = TGA + 1;
            end
        length(j) = (stop - start) + 3;
        start = 0;
    end
end

% Max, avg, and min protein lengths
avg_length = mean(length(length ~= 0));
max_length = max(length(length ~= 0));
min_length = min(length(length ~= 0));
proteins_segments = nnz(length);

% Print results
fprintf('\nTotal Protein-Coding Segments: %i\n', proteins_segments)
fprintf('Average Length: %.2f\n', avg_length)
fprintf('Maximum Length: %i\n', max_length)
fprintf('Minimum Length: %i\n\n', min_length)

%%=========================================================================
end
