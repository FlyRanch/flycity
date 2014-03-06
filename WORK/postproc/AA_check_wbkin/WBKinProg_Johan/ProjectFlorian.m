function ProjectFlorian( settings, pathDB )

    % Calculate a matrix with a column which contains the side-slip angle
    % in degrees and the yaw angle of the velocity vector in the global
    % reference frame.
    
    for i=1:size(pathDB.x,2)
        
       close all;

       start = find(isnan(pathDB.x(:,i))==0, 1 );
       stop = find(isnan(pathDB.x(:,i))==0, 1, 'last' );
    
        % Read pathDB
    
        t = pathDB.t(start:stop,1);
        
        u_f = pathDB.u_filt(start:stop,i);
        v_f = pathDB.v_filt(start:stop,i);
        
               
        b_beta = pathDB.b_beta(start:stop,i);
              
        % Calculate yaw angle of velocity vector
        
        V_yaw = real(atan2(v_f,u_f));
        
        figure()
        plot(t,radtodeg(b_beta))
        
        figure()
        plot(t,radtodeg(V_yaw))
        
        pause
        
            
    end


end

