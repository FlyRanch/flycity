function [ maneuver ] = select_maneuvers2(stroke_var,w_tresh,a_tresh)


    % Program that selects subsequent "maneuvering wingbeats":
    
    
    w_x = stroke_var.Omega_strk(1,:);
    
    w_y = stroke_var.Omega_strk(2,:);
    
    w_z = stroke_var.Omega_strk(3,:);
    
    
    a_x = stroke_var.A_strk(1,:);
    
    a_y = stroke_var.A_strk(2,:);
    
    a_z = stroke_var.A_strk(3,:);
    
    nr_wb = length(w_x);
    
    maneuver = {};
    
    maneuver.yaw_turns = zeros(1,nr_wb);
    maneuver.roll_turns = zeros(1,nr_wb);
    maneuver.pitch_turns = zeros(1,nr_wb);
    maneuver.a_x = zeros(1,nr_wb);
    maneuver.a_y = zeros(1,nr_wb);
    maneuver.a_z = zeros(1,nr_wb);

    
    for i = 2:(nr_wb-1)
        
        if w_x(i) > w_tresh
            maneuver.roll_turns(i) = i;
        end
        
        if w_y(i) > w_tresh
            maneuver.pitch_turns(i) = i;
        end
        
        if w_z(i) > w_tresh
            maneuver.yaw_turns(i) = i;
        end
        
        if w_x(i) < -w_tresh
            maneuver.roll_turns(i) = i;
        end
        
        if w_y(i) < -w_tresh
            maneuver.pitch_turns(i) = i;
        end
        
        if w_z(i) < -w_tresh
            maneuver.yaw_turns(i) = i;
        end
        
        if a_x(i) > a_tresh
            maneuver.a_x(i) = i;
        end
        
        if a_x(i) < -a_tresh
            maneuver.a_x(i) = i;
        end
        
        if a_y(i) > a_tresh
            maneuver.a_y(i) = i;
        end
        
        if a_y(i) < -a_tresh
            maneuver.a_y(i) = i;
        end
        
        if a_z(i) > a_tresh
            maneuver.a_z(i) = i;
        end
        
        if a_z(i) < -a_tresh
            maneuver.a_z(i) = i;
        end
        
    end
    
    for i = 2:(nr_wb-1)
        
        if maneuver.roll_turns(i) > 0
            
            if maneuver.roll_turns(i-1)==0 && maneuver.roll_turns(i+1)==0
                
                maneuver.roll_turns(i) = 0;
                
            elseif maneuver.roll_turns(i-1)==0 && maneuver.roll_turns(i+1) > 0
                
                if i > 2
                
                maneuver.roll_turns(i-1) = i-1;
                
                end
                
            end
            
        end
        
        if maneuver.pitch_turns(i) > 0
            
            if maneuver.pitch_turns(i-1)==0 && maneuver.pitch_turns(i+1)==0
                
                maneuver.pitch_turns(i) = 0;
                
            elseif maneuver.pitch_turns(i-1)==0 && maneuver.pitch_turns(i+1) > 0
                
                if i > 2
                
                maneuver.pitch_turns(i-1) = i-1;
                
                end
                
            end
            
        end
        
        if maneuver.yaw_turns(i) > 0
            
            if maneuver.yaw_turns(i-1)==0 && maneuver.yaw_turns(i+1)==0
                
                maneuver.yaw_turns(i) = 0;
                
            elseif maneuver.yaw_turns(i-1)==0 && maneuver.yaw_turns(i+1) > 0
                
                if i > 2
                
                maneuver.yaw_turns(i-1) = i-1;
                
                end
                
            end
            
        end
        
        if maneuver.a_x(i) > 0
            
            if maneuver.a_x(i-1)==0 && maneuver.a_x(i+1)==0
                
                maneuver.a_x(i) = 0;
                
            elseif maneuver.a_x(i-1)==0 && maneuver.a_x(i+1) > 0
                
                if i > 2
                
                maneuver.a_x(i-1) = i-1;
                
                end
                
            end
            
        end
        
        if maneuver.a_y(i) > 0
            
            if maneuver.a_y(i-1)==0 && maneuver.a_y(i+1)==0
                
                maneuver.a_y(i) = 0;
                
            elseif maneuver.a_y(i-1)==0 && maneuver.a_y(i+1) > 0
                
                if i > 2
                
                maneuver.a_y(i-1) = i-1;
                
                end
                
            end
            
        end
        
        if maneuver.a_z(i) > 0
            
            if maneuver.a_z(i-1)==0 && maneuver.a_z(i+1)==0
                
                maneuver.a_z(i) = 0;
                
            elseif maneuver.a_z(i-1)==0 && maneuver.a_z(i+1) > 0
                
                if i > 2
                
                maneuver.a_z(i-1) = i-1;
                
                end
                
            end
            
        end
        
    end

    
    for i = (nr_wb-1):-1:2
        
        if maneuver.roll_turns(i) > 0
            
            if maneuver.roll_turns(i+1)==0 && maneuver.roll_turns(i-1)==0
                
                maneuver.roll_turns(i) = 0;
                
            elseif maneuver.roll_turns(i+1)==0 && maneuver.roll_turns(i-1) > 0
                
                if i < (nr_wb-1)
                
                maneuver.roll_turns(i+1) = i+1;
                
                end
                
            end
            
        end
        
        if maneuver.pitch_turns(i) > 0
            
            if maneuver.pitch_turns(i+1)==0 && maneuver.pitch_turns(i-1)==0
                
                maneuver.pitch_turns(i) = 0;
                
            elseif maneuver.pitch_turns(i+1)==0 && maneuver.pitch_turns(i-1) > 0
                
                if i < (nr_wb-1)
                
                maneuver.pitch_turns(i+1) = i+1;
                
                end
                
            end
            
        end
        
        if maneuver.yaw_turns(i) > 0
            
            if maneuver.yaw_turns(i+1)==0 && maneuver.yaw_turns(i-1)==0
                
                maneuver.yaw_turns(i) = 0;
                
            elseif maneuver.pitch_turns(i+1)==0 && maneuver.yaw_turns(i-1) > 0
                
                if i < (nr_wb-1)
                
                maneuver.yaw_turns(i+1) = i+1;
                
                end
                
            end
            
        end
        
        if maneuver.a_x(i) > 0
            
            if maneuver.a_x(i+1)==0 && maneuver.a_x(i-1)==0
                
                maneuver.a_x(i) = 0;
                
            elseif maneuver.a_x(i+1)==0 && maneuver.a_x(i-1) > 0
                
                if i < (nr_wb-1)
                
                maneuver.a_x(i+1) = i+1;
                
                end
                
            end
            
        end
        
        if maneuver.a_y(i) > 0
            
            if maneuver.a_y(i+1)==0 && maneuver.a_y(i-1)==0
                
                maneuver.a_y(i) = 0;
                
            elseif maneuver.a_y(i+1)==0 && maneuver.a_y(i-1) > 0
                
                if i < (nr_wb-1)
                
                maneuver.a_y(i+1) = i+1;
                
                end
                
            end
            
        end
        
        if maneuver.a_z(i) > 0
            
            if maneuver.a_z(i+1)==0 && maneuver.a_z(i-1)==0
                
                maneuver.a_z(i) = 0;
                
            elseif maneuver.a_z(i+1)==0 && maneuver.a_z(i-1) > 0
                
                if i < (nr_wb-1)
                
                maneuver.a_z(i+1) = i+1;
                
                end
                
            end
            
        end
        
    end

end

