%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Homework 02: Main Matlab Script
%   [1] Neighbor Identification
%   [2] The Three Species Problem
%
%   Author: Parminder Singh
%   Date:   07/10/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clear Cache
clear all
close all
clc

%% Homework 01 Problems List
fprintf('Homework 01 Problems:\n\n')
fprintf('1: Neighbor Identification\n')
fprintf('2: The Three Species Problem\n\n')

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
M = input('\nEnter number of rows (M): ');
    while (M < 2) || abs(floor(M) - M) > 1e-10 %error check for M
        fprintf('ERROR: Value for M is not an integer greater than or equal to 2.\n');
        M = input('Enter M: ');
    end

N = input('Enter number of columns (N): ');
    while (N < 2) || abs(floor(N) - N) > 1e-10 %error check for N
        fprintf('ERROR: Value for N is not an integer greater than or equal to 2.\n');
        N = input('Enter N: ');
    end

MN = M * N; 
P = input('Enter Cell ID (P): '); %error check for P
    while P < 1 || P > (MN) || abs(floor(P) - P) > 1e-10 
        fprintf('ERROR: Value for P is not an integer within the bounds [1, M*N]\n');
        P = input('Enter P: ');
    end 


%left wall
if P > 1 && P < M  
    fprintf('\nLeft wall\n'); 
    fprintf('Cell ID: %d\n', P);
    fprintf('Neighbors: %d %d %d %d %d\n', P - 1 , P + 1, ((P + M) - 1), P + M, (P + M) + 1);

%right wall
elseif P <= (MN - 1) && P >= (M*(N - 1) + 1) +1 && P ~= MN && P ~= (M*(N - 1) - 1)
    fprintf('\nRight wall\n'); 
    fprintf('Cell ID: %d\n', P);
    fprintf('Neighbors: %d %d %d %d %d\n', (P - M) - 1, P - M, (P - M) + 1, P - 1, P + 1);
   
%bottom wall
elseif mod(P, M) == 0 && P ~= M && P ~= (M*N)
    fprintf('\nBottom wall\n'); 
    fprintf('Cell ID: %d\n', P);
    fprintf('Neighbors: %d %d %d %d %d\n', (P - M) - 1, P - M, P - 1, (P +M) - 1, P + M);
  
%top wall
elseif mod(P, (M - 1)) == 1 && P ~= 1 && P ~= ((M*(N-1)) + 1)
    fprintf('\nTop wall\n'); 
    fprintf('Cell ID: %d\n', P);
    fprintf('Neighbors: %d %d %d %d %d\n', P - M, (P - M) + 1, P + 1, P + M, (P + M) + 1);

%bottom left corner 
elseif P == M
    fprintf('\nBottom left corner\n'); 
    fprintf('Cell ID: %d\n', P);
    fprintf('Neighbors: %d %d %d\n', P - 1, (P + M) - 1, P + M);
    

%bottom right corner
elseif P == MN
    fprintf('\nBottom right corner\n'); 
    fprintf('Cell ID: %d\n', P);
    fprintf('Neighbors: %d %d %d\n', (P - M) - 1, P - M, P - 1);

%top left corner 
elseif P == 1
    fprintf('\nTop left corner\n'); 
    fprintf('Cell ID: %d\n', P);
    fprintf('Neighbors: %d %d %d\n', P + 1, P + M, (P + M) + 1);
    
%top right corner 
elseif P == ((M*(N-1)) + 1)
    fprintf('\nTop right corner\n'); 
    fprintf('Cell ID: %d\n', P);
    fprintf('Neighbors: %d %d %d\n', P - M, (P - M) + 1, P + 1);
    
%Interior
else
    fprintf('\nInterior\n'); 
    fprintf('Cell ID: %d\n', P);
    fprintf('Neighbors: %d %d %d %d %d %d %d %d\n', (P - M) - 1, P - M, (P - M) + 1, P -1, P +1, (P + M) - 1, P + M, (P + M) + 1);
end

%% Problem Number 2: The Three Speicies Problem ===========================
case 2

%Initial values
X_0 = 2;
Y_0 = 2.49;
Z_0 = 1.5;
T_final = 12;
delta_t = 0.005;

%print table with original values only
fprintf('Time\t X\t Y\t Z\n');
fprintf('0.00\t 2.00\t 2.49\t 1.50\n');

%calculate t_step
t_step = ceil(T_final/delta_t);

tic %use tic before for loop

for k = 1:t_step
    X_1 = X_0 + (delta_t*(0.75*X_0*(1 - (X_0/20)) - (1.5*X_0*Y_0) - (0.5*X_0*Z_0)));
    Y_1 = Y_0 + (delta_t*(Y_0*(1 - (Y_0/25)) - (0.75*X_0*Y_0) - (1.25*Y_0*Z_0)));
    Z_1 = Z_0 + (delta_t*(1.5*Z_0*(1 - (Z_0/30)) - (X_0*Z_0) - (Y_0*Z_0)));
    X_0 = X_1;
    Y_0 = Y_1;
    Z_0 = Z_1;
    t = k * delta_t;
    %Print rest of the table
    if mod(t, 0.5) == 0
        fprintf('%.2f\t %.2f\t %.2f\t %.2f\n', t, X_0, Y_0, Z_0);
    end
end

%time elapsed
tictoc = toc;
fprintf('Time elapsed: %f seconds\n', tictoc);

%%=========================================================================
end
