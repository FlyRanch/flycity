function [ theta_sd ] = Euler_sd(theta)

    % Calculates the mean without including a possible singularity in the
    % Euler angles:
        
    theta_plus = zeros(length(theta),1);
    
    % Calculates the mean without including a possible singularity in the
    % Euler angles:
    
    theta_min = zeros(length(theta),1);
    
    theta_plus = zeros(length(theta),1);
 
    N = length(theta);
    
    for i = 1:length(theta)
        
        if theta(i) > 0
            
            theta_plus(i) = theta(i);

            
        elseif theta(i) <0
            
            theta_plus(i) = 2*pi+theta(i);

            
        end
        
    end
    
    theta_mean = sum(theta_plus)/N;
    
    % Calculate the standard deviation without a possible singularity due
    % to the Euler angles:
    
    N = length(theta);
    
    theta_sd = sqrt((1/(N-1)).*sum((abs(theta)-abs(theta_mean)).^2));



end
