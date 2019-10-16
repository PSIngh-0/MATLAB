%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Final Project: Spinodal Decomposition
% 
% Author: Parminder Singh
% Date:   08/14/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clear Cache
clear all
close all
clc

%% Main Script
fprintf('Final Problems:\n\n')
fprintf('1: Five-Point Stencil with First Order Runge-Kutta Approximation\n')
fprintf('2: Nine-Point Stencil with First Order Runge-Kutta Approximation\n')
fprintf('3: Five-Point Stencil with Second Order Runge-Kutta Approximation\n\n')

problem = input('Please enter the problem number1: \n');
    while (problem < 1 || problem > 3) %error check for problem input
        fprintf('\nERROR: Invalid problem number! Please try again.\n\n');
        fprintf('1: Five-Point Stencil with First Order Runge-Kutta Approximation\n')
        fprintf('2: Nine-Point Stencil with First Order Runge-Kutta Approximation\n')
        fprintf('3: Five-Point Stencil with Second Order Runge-Kutta Approximation\n\n')
        problem = input('Please enter the problem number: \n');
    end
    
% Define given constant variables
h = 1;
a = 1;
b = 1;
gamma = 1;
D = 3;

row = 150;
col = 100;

% Time information split up and named using given values from parts 1, 2,
% 3, etc...
t_0 = 0;
t_f = 10;
t_s = 20;
dt_1 = 1e-3;
dt_2 = 1e-4;
dt_3 = 1e-3;
nt_1 = (t_f - t_0)/dt_1;
nt_2 = (t_f - t_0)/dt_2;
nt_3 = (t_s - t_0)/dt_3;
frames = 301;

rng('shuffle');
rng('default')
phi = rand(row,col);
for i = 1:1:row
    for j = 1:1:col
        p = rand; 
         if phi(i, j) < 0.5
            phi(i, j) = -1;
           
        else
            phi(i, j) = 1;
        end
    end
end




switch(problem)

    % First Order Runge Kutta Method with five-point stencil
    case 1
        
        vidfile = VideoWriter('Project_Problem_1.mp4','MPEG-4');
        vidfile.FrameRate = 30;
        open(vidfile);
        
        figure(1)
        imagesc(phi)
        colorbar
        caxis([-1, 1])
        title('Spinodal Decomposition')
        drawnow
        
        for k = 1:1:nt_1
            v_5_1 = Laplacian_2D(phi, h, 5);
            u_1 = ((((b^4)*(phi.^3)) - (a*(b^2)*phi)) - gamma*v_5_1);
            v_5_intermediate_1 = Laplacian_2D(u_1, h, 5);
            d_phi_1 = D*v_5_intermediate_1;

            C1_1 = dt_1*(d_phi_1);
            d_phi_new_1 = phi + C1_1;

            phi = d_phi_new_1;
            
            t_intermediate = t_f/dt_1;
            
            t_inter = t_intermediate/frames;
            
            t_final = floor(k/t_inter);
            
            frame_video = 0;
            
            if t_final <= frames
                t_final = (t_f/frames)*t_final;
                imagesc(phi);
                colorbar
                caxis([-1, 1])
                title('Spinodal Decomposition')
                drawnow
                frame_video = frame_video + 1;
                Fr(frame_video) = getframe(gcf);
                writeVideo(vidfile, Fr);                
            end
            
        end 

    % First Order Runge Kutta Method with nine-point stencil    
    case 2
        vidfile = VideoWriter('Project_Problem_2.mp4','MPEG-4');
        vidfile.FrameRate = 30;
        open(vidfile);
        
        figure(1)
        imagesc(phi)
        colorbar
        caxis([-1, 1])
        title('Spinodal Decomposition')
        drawnow
        
        
        for k = 1:1:nt_2
            v_9_2 = Laplacian_2D(phi, h, 9);
            u_2 = ((((b^4)*(phi.^3)) - (a*(b^2)*phi)) - gamma*v_9_2);
            v_9_intermediate_2 = Laplacian_2D(u_2, h, 9);
            d_phi_2 = D*v_9_intermediate_2;

            C1_2 = dt_2*(d_phi_2);
            d_phi_new_2 = phi + C1_2;

            phi = d_phi_new_2;
            
            t_intermediate = t_f/dt_1;
            
            t_inter = t_intermediate/frames;
            
            t_final = floor(k/t_inter);
            
            frame_video = 0;
            
            if t_final <= frames
                t_final = (t_f/frames)*t_final;
                imagesc(phi);
                colorbar
                caxis([-1, 1])
                title('Spinodal Decomposition')
                drawnow
                frame_video = frame_video + 1;
                Fr(frame_video) = getframe(gcf);
                writeVideo(vidfile, Fr);                 
            end

        end 

    % Second Order Runge Kutta Method with nine-point stencil    
    case 3
        
        vidfile = VideoWriter('Project_Problem_3.mp4','MPEG-4');
        vidfile.FrameRate = 30;
        open(vidfile);        
        
        
        figure(1)
        imagesc(phi)
        colorbar
        caxis([-1, 1])
        title('Spinodal Decomposition t = 0')
        drawnow
        
        
        for k = 1:1:nt_3
            v_9_3 = Laplacian_2D(phi, h, 5);
            u_3 = ((((b^4)*(phi.^3)) - (a*(b^2)*phi)) - gamma*v_9_3);
            v_9_intermediate_3 = Laplacian_2D(u_3, h, 5);
            d_phi_3 = D*v_9_intermediate_3;

            C1_3 = dt_3*d_phi_3;
            d_phi_new_3 = phi + 0.5*C1_3;

            v_9_3_2 = Laplacian_2D(d_phi_new_3, h, 5);
            u_3_2 = ((((b^4)*(d_phi_new_3.^3)) - (a*(b^2)*d_phi_new_3)) - gamma*v_9_3_2);
            v_9_intermediate_3_2 = Laplacian_2D(u_3_2, h, 5);
            d_phi_3_2 = D*v_9_intermediate_3_2;            
            
            C2_3 = dt_3*(d_phi_3_2);
            d_phi_new_3_2 = phi + C2_3;         
            
            phi = d_phi_new_3_2;
            
            t_f_3 = dt_3*k;
            
            frame_video = 0;
            
            if mod(t_f_3, 5) == 0
                imagesc(d_phi_new_3_2);
                colorbar
                caxis([-1, 1])
                title(['Spinodal Decomposition t = ' num2str(t_f_3) '.'])
                drawnow
                frame_video = frame_video + 1;
                Fr(frame_video) = getframe(gcf);
                writeVideo(vidfile, Fr);  
            end
            

        end
    
end