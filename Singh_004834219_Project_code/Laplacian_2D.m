
function [ v ] = Laplacian_2D(phi, h, stencil)

row = 150;
col = 100;

v = zeros(row, col);

switch stencil
    case 5
        for i = 1:1:row 
            for j = 1:1:col
                phi_left = j - 1;
                phi_right = j + 1;
                phi_down = i + 1;
                phi_up = i - 1;
                if phi_left < 1
                    phi_left = col;
                end
                if phi_right > col
                    phi_right = 1;
                end
                if phi_down > row
                    phi_down = 1;
                end
                if phi_up < 1
                    phi_up = row;
                end
                
                approximation = phi(i, phi_left) + phi(i, phi_right) + ...
                                phi(phi_down, j) + phi(phi_up, j) - ...
                                (4*phi(i, j));
                
                v_intermediate = (1/h^2)*(approximation);
                v(i, j) = v_intermediate;
            end
        end 
        
                
                
    case 9
         for i = 1:1:row
            for j = 1:1:col
                phi_left = j - 1;
                phi_right = j + 1;
                phi_down = i + 1;
                phi_up = i - 1;
                if phi_left < 1
                    phi_left = col;
                end
                if phi_right > j
                    phi_right = 1;
                end
                if phi_down > i
                    phi_down = 1;
                end
                if phi_up < 1
                    phi_up = row;
                end
                
                approximation = phi(phi_up, phi_left) + phi(phi_down, phi_left) + ...
                                phi(phi_up, phi_right) + phi(phi_down, phi_right) +...
                                (4*phi(i, phi_left)) + (4*phi(i, phi_right)) + ...
                                (4*phi(phi_down, j)) + (4*phi(phi_up, j)) - ...
                                (20*phi(i, j));
                
                v_intermediate = (1/(6*h^2))*(approximation); 
                v(i, j) = v_intermediate;
            end 
         end 
      
end

end 