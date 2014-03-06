function Histogram_Floris(settings, pathDB)


    %Make histogram of sideslip angle vs absolute horizontal velocity
    
    M = size(pathDB.x,2);
    
    beta = nan(size(pathDB.x,1),M);
    
    vel = nan(size(pathDB.x,1),M);
    
    
    
    for i=1:M
        
        % start and stop point for the measurements
        start = find(isnan(pathDB.x(:,i))==0, 1 );
        stop = find(isnan(pathDB.x(:,i))==0, 1, 'last' );
        
        %Remove nan-values from necessary variables in pathDB:        
    
        u = pathDB.u_filt(start:stop,i);
        v = pathDB.v_filt(start:stop,i);
        
        yaw = pathDB.b_yaw(start:stop,i);
        
        % Calculate the horizontal velocity vector and the horizontal
        % heading vector:
        
        N = stop-start+1;
        
        vel_vec = zeros(2,N);
        
        head_vec = zeros(2,N);
        
        for j = 1:N
            
            vel_vec(:,j) = [u(j); v(j)];
            
            head_vec(:,j) = [cos(yaw(j)); sin(yaw(j))];

            beta(start-1+j,i) = radtodeg(real(acos(dot(vel_vec(:,j),head_vec(:,j))/(norm(vel_vec(:,j))*norm(head_vec(:,j))))));
            
            vel(start-1+j,i) = norm(vel_vec(:,j));
                        
        end
        
        clear start stop
        
        
    end
    
    
    start_t = find(isnan(beta(:,1)) == 0, 1 );
    stop_t = find(isnan(beta(:,1)) == 0, 1, 'last' );
    
    beta_hist_t = beta(start_t:stop_t,1);
    
    vel_hist_t =  vel(start_t:stop_t,1);


    
    % Make a histogram of all the sequences
    
    % List all the beta's and velocities for a histogram:
    
    for i = 2:M
        
        clear start_t stop_t
        
        start_t = find(isnan(beta(:,i)) == 0, 1 );
        %stop_t = find(isnan(beta(:,i)) == 0, 1, 'last' );
        stop_t = 2795;

        beta_hist = [ beta_hist_t; beta(start_t:stop_t,i)];
    
        vel_hist = [ vel_hist_t; vel(start_t:stop_t,i)];
        
        beta_hist_t = beta_hist;
        
        vel_hist_t = vel_hist;
        
        
    end

    figure()
    scatter(vel_hist,beta_hist)

     [N_histx,x_hist] = hist(vel_hist,15);

     [N_histy, y_hist] = hist(beta_hist,15);
    
    [pdf_z,B] = hist3([vel_hist beta_hist], [15 15]);
    
    figure()
    surf(x_hist,y_hist,pdf_z)
    title('3D histogram speed v.s. sideslip angle')
    xlabel('Horizontal velocity [mm/s]')
    ylabel('Sideslip angle \beta [deg]')
    zlabel('Nr of counts')

    
    
end

