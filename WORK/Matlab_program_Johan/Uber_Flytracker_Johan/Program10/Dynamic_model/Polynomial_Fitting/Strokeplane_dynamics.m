function [ stroke_var ] = Strokeplane_dynamics( settings, pathDB, seq_nr )


    % Return the following variables in the strokeplane reference frame:
    
    % Vn, Vt, wx, wy, wz
    
    % An, At, w_dot_x, w_dot_y, w_dot_z
    
    %----------------------------------------------------------------------
    
    stroke_var = {};
    
    start = find(isnan(pathDB.x(:,seq_nr))==0, 1 );
    stop = find(isnan(pathDB.x(:,seq_nr))==0, 1, 'last' );
    
    t = pathDB.t(start:stop);
    
    dt = pathDB.t(2)-pathDB.t(1);
    
    
    nr_wb = find(isnan(pathDB.wingbeat_time(:,1,seq_nr))==0, 1, 'last' );
    
    q_body = zeros(4,nr_wb);
    
    Omega_strk = zeros(3,nr_wb);
    
    Omega_dot_strk = zeros(3,nr_wb);
    
    V_strk = zeros(3,nr_wb);
    
    A_strk = zeros(3,nr_wb);
    
    
    for i = 1:nr_wb
        
        a = find(isnan(pathDB.wingbeat_time(i,:,seq_nr))==0, 1, 'last' );
        
        b1 = start-1+pathDB.wingbeat_time(i,1:(a-1),seq_nr);
        
        b2 = start-1+pathDB.wingbeat_time(i,2:a,seq_nr);
        
        q1_wb = pathDB.qb1_filt(start-1+pathDB.wingbeat_time(i,1:a,seq_nr),seq_nr);
        q2_wb = pathDB.qb2_filt(start-1+pathDB.wingbeat_time(i,1:a,seq_nr),seq_nr);
        q3_wb = pathDB.qb3_filt(start-1+pathDB.wingbeat_time(i,1:a,seq_nr),seq_nr);
        q4_wb = pathDB.qb4_filt(start-1+pathDB.wingbeat_time(i,1:a,seq_nr),seq_nr);
        
        q_wb_mean = q_avg(q1_wb, q2_wb ,q3_wb, q4_wb);
        
                
        u_body = pathDB.u_body_mean(i,1,seq_nr);
        
        v_body = pathDB.v_body_mean(i,1,seq_nr);
        
        w_body = pathDB.w_body_mean(i,1,seq_nr);
        
        wx_body = pathDB.omegax_body_mean(i,1,seq_nr);
        
        wy_body = pathDB.omegay_body_mean(i,1,seq_nr);
        
        wz_body = pathDB.omegaz_body_mean(i,1,seq_nr);
        
               
        ax_body = pathDB.ax_body_mean(i,1,seq_nr);
        
        ay_body = pathDB.ay_body_mean(i,1,seq_nr);
        
        az_body = pathDB.az_body_mean(i,1,seq_nr);
        
        w_dot_x_body = mean((pathDB.b_omega1(b2,seq_nr)-pathDB.b_omega1(b1,seq_nr))./dt);
        
        w_dot_y_body = mean((pathDB.b_omega2(b2,seq_nr)-pathDB.b_omega2(b1,seq_nr))./dt);
        
        w_dot_z_body = mean((pathDB.b_omega3(b2,seq_nr)-pathDB.b_omega3(b1,seq_nr))./dt);
        
        
        % Transfer the variables to the strokeplane
        
        beta = -(55/180)*pi;    
    
        R_beta = [cos(beta) 0 -sin(beta); ...
               0 1 0; ...
               sin(beta) 0 cos(beta)]; 
        
        q_body(:,i) = q_wb_mean; 
           
        V_strk(:,i) = R_beta*[u_body; v_body; w_body];
        
        Omega_strk(:,i) = R_beta*[wx_body; wy_body; wz_body];
        
        A_strk(:,i) = R_beta*[ax_body; ay_body; az_body];
        
        Omega_dot_strk(:,i) = R_beta*[w_dot_x_body; w_dot_y_body; w_dot_z_body];
                
        
    end

    stroke_var.q_body = q_body;
    
    stroke_var.Omega_strk = Omega_strk;
    
    stroke_var.Omega_dot_strk = Omega_dot_strk;
    
    stroke_var.V_strk = V_strk;
    
    stroke_var.A_strk = A_strk;

end

