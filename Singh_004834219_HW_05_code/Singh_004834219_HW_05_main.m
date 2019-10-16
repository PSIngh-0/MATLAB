%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Homework 05: Main Matlab Script
% [1] The Shared Birthday Problem
% [2] Random Walk Collisions
%
% Author: Parminder Singh
% Date:   07/29/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clear Cache
clear all
close all
clc

%% Homework 04 Problems List
fprintf('Homework 05 Problems:\n\n')
fprintf('1: The Shared Birthday Problem\n')
fprintf('2: Random Walk Collisions\n\n')

problem = input('Please enter the problem number: \n');
    while (problem <= 0 || problem > 2) %error check for problem input
        fprintf('\nERROR: Invalid problem number! Please try again.\n\n');
        fprintf('1: The Shared Birthday Problem\n')
        fprintf('2: Random Walk Collisions\n\n')
        problem = input('Please enter the problem number: \n');
    end
switch(problem)

    
%% Problem Number 1: The Shared Birthday Problem ==========================
case 1
    
rng('shuffle');    
trials = 1e4;
birthdays = zeros(trials, 1);

% For loop which initiates variables used for the while loop within this
% loop and also finds a random integer between 1 and 365
for k = 1:1:trials
    
    x = 1;
    birthday(x) = randi([1 365]);
    
    j = 0;

    
    % While loop which creates random integers between 1 and 365 and checks
    % if they are within 7 1days
    while j == 0
        x = x + 1;        
        birthday(x) = randi([1 365]); 
        
        % For loop which checks if the birthday group is within 7 days and
        % stop it if it is.
        for i = 1:1:(x - 1)
            
            if  abs(birthday(x) - birthday(i)) < 8 || abs(birthday(x) - birthday(i)) > 358
                j = 1; 
                birthdays(k) = x;
            end 
            
        end 
        
    end
    
end

median = median(birthdays);
fprintf('\nMedian Number of People = %i \n\n', median)

histogram(birthdays)
title('Distribution of groups of people who have their birthdays within a week')
xlabel('Number of people in a group')
ylabel('Number of trials with the same amount of people in a group')

%% Problem Number 2: Random Walk Collisions ===============================
case 2

rng('shuffle');

trials = 5e3;
boundary = [5, -5, 5, -5];

% Set initial postions for particles 1 and 2
x1 = -5;
y1 = 0;
x2 = 5;
y2 = 0;

moveCount = 0;
maxCount = 5e3;
actualmoveCount = zeros(trials, 1);

% For loop responsible for using move function to randomly assign positions
% to x1, x2, y1, and y2
for k = 1:1:trials
    
    % While loop which stops once collision conditions are met or the
    % amount of moves exceeds the maximum number of allowed moves
    while (moveCount < maxCount) && (x1 ~= x2 || y1 ~= y2)  
        [x1, y1] = rand_walk(x1, y1, boundary);
        [x2, y2] = rand_walk(x2, y2, boundary);
        moveCount = moveCount + 1;
       
        while k == 1 && (moveCount < maxCount) && (x1 ~= x2 || y1 ~= y2)
            % Values for initial points
            x = -5;
            y = 0;
            x0 = 5;
            y0 = 0;

            % Generate random step for both A and B
            [x1, y1] = rand_walk(x1, y1, boundary);

            [x2, y2] = rand_walk(x2, y2, boundary);

            % Initial and Final squares for A
            xk_val = [x - 0.5, x + 0.5, x + 0.5, x - 0.5];
            yk_val = [y - 0.5, y - 0.5, y + 0.5, y + 0.5];

            xkp_val = [x1 - 0.5, x1 + 0.5, x1 + 0.5, x1 - 0.5];
            ykp_val = [y1 - 0.5, y1 - 0.5, y1 + 0.5, y1 + 0.5];

            % Initial and Final squares for B
            x2k_val = [x0 - 0.5, x0 + 0.5, x0 + 0.5, x0 - 0.5];
            y2k_val = [y0 - 0.5, y0 - 0.5, y0 + 0.5, y0 + 0.5];

            x2kp_val = [x2 - 0.5, x2 + 0.5, x2 + 0.5, x2 - 0.5];
            y2kp_val = [y2 - 0.5, y2 - 0.5, y2 + 0.5, y2 + 0.5];

            % Plot with both A and B
            figure(1)

            hold on

            set(gca, 'xtick', -5.5:1:5.5)
            set(gca, 'ytick', -5.5:1:5.5)

            grid on
            xlim([-5.5 5.5])
            ylim([-5.5 5.5])

            axis square
            fill(xk_val, yk_val, 'r') % Intial for A
                
            fill(xkp_val, ykp_val, 'b') % Final for A
                
            fill(x2k_val, y2k_val, 'r') % Intial for B
                
            fill(x2kp_val, y2kp_val, 'b') % Final for A
            title('2D Random Walk')

            hold off
        end
    end 
    
    actualmoveCount(k) = moveCount;
    
    % Define initial positions and move counter again for the loop to run 
    % again with initial values
    x1 = -5;
    y1 = 0;
    x2 = 5;
    y2 = 0;
    moveCount = 0;
    
end 

% Calculate median value of moves for the trials
x = median(actualmoveCount);
median = ceil(x);
fprintf('Median = %i\n', median)


%%=========================================================================
end
