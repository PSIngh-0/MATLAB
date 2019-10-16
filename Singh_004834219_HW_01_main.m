%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Homework 01: Main Matlab Script
%   [1] Oblate Spherical Calculations
%   [2] Ellipse Perimeter Calculation
%
%   Author: Parminder Singh
%   Date:   06/28/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clear Cache
clear all
close all
clc

%% Homework 01 Problems List
fprintf('Homework 01 Problems:\n\n')
fprintf('1: Oblate Spherical Calculations\n')
fprintf('2: Ellipse Perimeter Calculation\n\n')

problem = input('Please enter the problem number: \n');

switch(problem)
    
%% Problem Number 1: Oblate Spherical Calculations=========================
    
case 1
% Input values for r_1 and r_2
r_1 = input('Enter r_1: ');
r_2 = input('Enter r_2: ');

% Given equations
gamma = acos((r_2)/(r_1));
A_actual = 2*pi*((r_1)^2 + (((r_2)^2)/sin(gamma))*(log((cos(gamma))/(1 - sin(gamma)))));
A_approximate = 4*pi*((r_1 + r_2)/2)^2;
%percent error equation
percent_error = ((A_actual - A_approximate)/A_actual)*100; 

% Solving for inputted r_1 and r_2 values and printing results along with
% error
if (r_2 > r_1) || (r_1 < 0) || (r_2 < 0)
    
    fprintf('Please rerun the script and make sure to enter a value of r_2 less than r_1 and make sure to only enter real values greater than zero for r_1 and r_2\n') 
   
else
    fprintf('The actual surface area with inputted data is: %10f\n', A_actual);
    fprintf('The approximated surface area with inputted data is: %10f\n', A_approximate);
    fprintf('The percent error between the two methods of calculation is: %10f\n\n', percent_error);
    
end    

% Earth data along with error
r_e_1 = 6378.137;
r_e_2 = 6356.752;
gamma_e = acos((r_e_2)/(r_e_1));
A_e_actual = 2*pi*((r_e_1)^2 + (((r_e_2)^2)/sin(gamma_e))*(log((cos(gamma_e))/(1 - sin(gamma_e)))));
A_e_approximate = 4*pi*((r_e_1 + r_e_2)/2)^2;
percent_error_e = ((A_e_actual - A_e_approximate)/A_e_actual)*100;

% Printing of results using Earth data
fprintf('The actual surface area with Earth data is: %10f\n', A_e_actual);
fprintf('The approximated surface area with Earth data is: %10f\n', A_e_approximate);
fprintf('The percent error between the two methods of calculation is: %10f\n', percent_error_e);      
%==========================================================================

%% Problem Number 2: Ellipse Perimeter Calculation=========================

case 2
% Input values for a and b
a = input('Enter a: ');
b = input('Enter b: ');

% Perimeter calculations
h = ((a-b)/(a+b))^2;
P_1 = pi*(a+b);
P_2 = pi*sqrt(2*((a^2)+(b^2)));
P_3 = pi*sqrt(2*((a^2)+(b^2))-((a-b)^2)/2);
P_4 = pi*(a+b)*((1+(h/8))^2);
P_5 = pi*(a+b)*(1+((3*h)/(10+sqrt(4-3*h))));
P_6 = pi*(a+b)*(64-3*h^2)/(64-16*h);
P_7 = pi*(a+b)*(256-48*h-21*(h^2))/(256-112*h+3*(h^2));
P_8 = pi*(a+b)*((3-sqrt(1-h))/2);

%Perimeters calulcated with the given formulas
if (a < 0 || b < 0)
    
    fprintf('Please rerun the script and enter real values greater than zero for a and b\n') 
   
else
    fprintf('\nPerimeters calculated with the given formulas:\n')
    fprintf('P_1 = %f\n',P_1);
    fprintf('P_2 = %f\n',P_2);
    fprintf('P_3 = %f\n',P_3);
    fprintf('P_4 = %f\n',P_4);
    fprintf('P_5 = %f\n',P_5);
    fprintf('P_6 = %f\n',P_6);
    fprintf('P_7 = %f\n',P_7);
    fprintf('P_8 = %f\n',P_8);
    fprintf('h = %f\n',h);
end
%==========================================================================
end
    