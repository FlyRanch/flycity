function Wingkinematics_Legendre(settings, pathDB)


    % Use a Legendre polynomial fit to describe the wingkinematics in a
    % minimal number of parameters:
    
    %     n_pol_theta = 13; % Order of used polynomials
%     n_pol_eta = 16; % Order of used polynomials
%     n_pol_phi = 8; % Order of used polynomials
%     
    n_pol_theta = 12; % Order of used polynomials
    n_pol_eta = 14; % Order of used polynomials
    n_pol_phi = 10; % Order of used polynomials
%     n_pol_theta = 9; % Order of used polynomials
%     n_pol_eta = 9; % Order of used polynomials
%     n_pol_phi = 7; % Order of used polynomials


    n_pol = {};
    
    n_pol.theta = n_pol_theta;
    
    n_pol.eta = n_pol_eta;
    
    n_pol.phi = n_pol_phi;
    

    
    dt = pathDB.t(2)-pathDB.t(1);
    
    nr_of_seq = size(pathDB.x,2);
    
    Maneuver_Matrix = zeros(1,13);
    
    a_sym_maneuver  = {};
    
    a_dev_maneuver = {};
    
    
    a_sym_maneuver.theta_L = zeros(1,n_pol_theta+1);
    
    a_sym_maneuver.theta_R = zeros(1,n_pol_theta+1);
    
    a_sym_maneuver.eta_L = zeros(1,n_pol_eta+1);
    
    a_sym_maneuver.eta_R = zeros(1,n_pol_eta+1);
    
    a_sym_maneuver.phi_L = zeros(1,n_pol_phi+1);
    
    a_sym_maneuver.phi_R = zeros(1,n_pol_phi+1);
    
    
    
    a_dev_maneuver.theta_L = zeros(1,n_pol_theta+1);
    
    a_dev_maneuver.theta_R = zeros(1,n_pol_theta+1);
    
    a_dev_maneuver.eta_L = zeros(1,n_pol_eta+1);
    
    a_dev_maneuver.eta_R = zeros(1,n_pol_eta+1);
    
    a_dev_maneuver.phi_L = zeros(1,n_pol_phi+1);
    
    a_dev_maneuver.phi_R = zeros(1,n_pol_phi+1);
    
    
    A_theta_L = 0;
    A_eta_L = 0;
    A_phi_L = 0;
    
    B_theta_L = 0;
    B_eta_L = 0;
    B_phi_L = 0;
    
    C_theta_L = 0;
    C_eta_L = 0;
    C_phi_L = 0;
    
    D_theta_L = 0;
    D_eta_L = 0;
    D_phi_L = 0;
    
    
    A_theta_R = 0;
    A_eta_R = 0;
    A_phi_R = 0;
    
    B_theta_R = 0;
    B_eta_R = 0;
    B_phi_R = 0;
    
    C_theta_R = 0;
    C_eta_R = 0;
    C_phi_R = 0;
    
    D_theta_R = 0;
    D_eta_R = 0;
    D_phi_R = 0;
    
    
    for i =1:nr_of_seq
    
    seq_nr = i
    
    start = find(isnan(pathDB.x(:,seq_nr))==0, 1 );
    stop = find(isnan(pathDB.x(:,seq_nr))==0, 1, 'last' );
    
    nr_wb = find(isnan(pathDB.wingbeat_time(:,1,seq_nr))==0, 1, 'last' );
    
    clear a_fit a_avg


    [a_fit,a_avg,f_avg,down_up_avg,trigger_wb,down_up_ratio] = Standard_wingbeat( settings, pathDB, i, n_pol_theta, n_pol_eta, n_pol_phi);

   
    
%     figure()
%     hold on
%     for k = 1:nr_wb
%     [ t_k, X_k ] = Wingbeat_Legendre_matrix( n_pol_eta, down_up_ratio(k), 200, 0, 1, 0 );
%     plot(t_k,X_k*[a_fit.eta_L1(:,k); a_fit.eta_L2(:,k)],'Color',[0.5 0.5 0.5])
%     end
%     [ t_k, X_k ] = Wingbeat_Legendre_matrix( n_pol_eta, down_up_avg, 200, 0, 1, 0 );
%     plot(t_k,X_k*[a_avg.eta_L1; a_avg.eta_L2],'r')
%     hold off
%     
%     
%     figure()
%     hold on
%     for k = 1:nr_wb
%     [ t_k, X_k ] = Wingbeat_Legendre_matrix( n_pol_eta, down_up_ratio(k), 200, 0, 1, 0 );
%     plot(t_k,X_k*[a_fit.eta_R1(:,k); a_fit.eta_R2(:,k)],'Color',[0.5 0.5 0.5])
%     end
%     [ t_k, X_k ] = Wingbeat_Legendre_matrix( n_pol_eta, down_up_avg, 200, 0, 1, 0 );
%     plot(t_k,X_k*[a_avg.eta_R1; a_avg.eta_R2],'r')
%     hold off
% 
%     figure()
%     hold on
%     for k = 1:nr_wb
%     [ t_k, X_k ] = Wingbeat_Legendre_matrix( n_pol_eta, down_up_ratio(k), 200, 0, 1, 1 );
%     plot(t_k,X_k*[a_fit.eta_L1(:,k); a_fit.eta_L2(:,k)],'Color',[0.5 0.5 0.5])
%     end
%     [ t_k, X_k ] = Wingbeat_Legendre_matrix( n_pol_eta, down_up_avg, 200, 0, 1, 1 );
%     plot(t_k,X_k*[a_avg.eta_L1; a_avg.eta_L2],'r')
%     hold off
%     
%     
%     figure()
%     hold on
%     for k = 1:nr_wb
%     [ t_k, X_k ] = Wingbeat_Legendre_matrix( n_pol_eta, down_up_ratio(k), 200, 0, 1, 1 );
%     plot(t_k,X_k*[a_fit.eta_R1(:,k); a_fit.eta_R2(:,k)],'Color',[0.5 0.5 0.5])
%     end
%     [ t_k, X_k ] = Wingbeat_Legendre_matrix( n_pol_eta, down_up_avg, 200, 0, 1, 1 );
%     plot(t_k,X_k*[a_avg.eta_R1; a_avg.eta_R2],'r')
%     hold off
%     
%     pause

    
    if trigger_wb >= 3 && trigger_wb <(nr_wb-1)


    [a_sym, a_dev] = Maneuver_wingbeat(a_fit,a_avg,n_pol_theta,n_pol_eta,n_pol_phi,down_up_ratio,down_up_avg);
    
    
    
    [ stroke_var ] = Strokeplane_dynamics( settings, pathDB, i );

    
    
    dev_theta_L = zeros(100*nr_wb,1);
    
    dev_theta_R = zeros(100*nr_wb,1);
    
    dev_eta_L = zeros(100*nr_wb,1);
    
    dev_eta_R = zeros(100*nr_wb,1);
    
    dev_phi_L = zeros(100*nr_wb,1);
    
    dev_phi_R = zeros(100*nr_wb,1);
    
    dev_theta_dot_L = zeros(100*nr_wb,1);
    
    dev_theta_dot_R = zeros(100*nr_wb,1);
    
    dev_eta_dot_L = zeros(100*nr_wb,1);
    
    dev_eta_dot_R = zeros(100*nr_wb,1);
    
    dev_phi_dot_L = zeros(100*nr_wb,1);
    
    dev_phi_dot_R = zeros(100*nr_wb,1);
    
    time_loc = pathDB.wingbeat_time(1:nr_wb,1,seq_nr);
    
    time = pathDB.t(start-1+time_loc);
    
    time2 = zeros(100*nr_wb,1);
    
  
    for k = 1:nr_wb
        
        if k < nr_wb
    
            [ ~ , X_theta ] = Wingbeat_Legendre_matrix( n_pol_theta, down_up_ratio(k), 100, time(k), time(k+1), 0 );

            [ ~ , X_eta ] = Wingbeat_Legendre_matrix( n_pol_eta, down_up_ratio(k), 100, time(k), time(k+1), 0 );

            [ t_temp , X_phi ] = Wingbeat_Legendre_matrix( n_pol_phi, down_up_ratio(k), 100, time(k), time(k+1), 0 );

            time2(((k-1)*100+1):(k*100)) = t_temp;

            [ ~ , X_dot_theta ] = Wingbeat_Legendre_matrix( n_pol_theta, down_up_ratio(k), 100, time(k), time(k+1), 1 );

            [ ~ , X_dot_eta ] = Wingbeat_Legendre_matrix( n_pol_eta, down_up_ratio(k), 100, time(k), time(k+1), 1 );

            [ ~ , X_dot_phi ] = Wingbeat_Legendre_matrix( n_pol_phi, down_up_ratio(k), 100, time(k), time(k+1), 1 );

            dev_theta_L(((k-1)*100+1):(k*100)) = X_theta*[a_dev.theta_L1(:,k); a_dev.theta_L2(:,k)];

            dev_theta_R(((k-1)*100+1):(k*100)) = X_theta*[a_dev.theta_R1(:,k); a_dev.theta_R2(:,k)];

            dev_eta_L(((k-1)*100+1):(k*100)) = X_eta*[a_dev.eta_L1(:,k); a_dev.eta_L2(:,k)];

            dev_eta_R(((k-1)*100+1):(k*100)) = X_eta*[a_dev.eta_R1(:,k); a_dev.eta_R2(:,k)];

            dev_phi_L(((k-1)*100+1):(k*100)) = X_phi*[a_dev.phi_L1(:,k); a_dev.phi_L2(:,k)];

            dev_phi_R(((k-1)*100+1):(k*100)) = X_phi*[a_dev.phi_R1(:,k); a_dev.phi_R2(:,k)];

            dev_theta_dot_L(((k-1)*100+1):(k*100)) = X_dot_theta*[a_dev.theta_L1(:,k); a_dev.theta_L2(:,k)];

            dev_theta_dot_R(((k-1)*100+1):(k*100)) = X_dot_theta*[a_dev.theta_R1(:,k); a_dev.theta_R2(:,k)];

            dev_eta_dot_L(((k-1)*100+1):(k*100)) = X_dot_eta*[a_dev.eta_L1(:,k); a_dev.eta_L2(:,k)];

            dev_eta_dot_R(((k-1)*100+1):(k*100)) = X_dot_eta*[a_dev.eta_R1(:,k); a_dev.eta_R2(:,k)];

            dev_phi_dot_L(((k-1)*100+1):(k*100)) = X_dot_phi*[a_dev.phi_L1(:,k); a_dev.phi_L2(:,k)];

            dev_phi_dot_R(((k-1)*100+1):(k*100)) = X_dot_phi*[a_dev.phi_R1(:,k); a_dev.phi_R2(:,k)];
        
        elseif k == nr_wb
            
            t_end = find(isnan(pathDB.wingbeat_time(nr_wb,:,seq_nr))==0, 1, 'last' );
            
            [ ~ , X_theta ] = Wingbeat_Legendre_matrix( n_pol_theta, down_up_ratio(k), 100, time(k), time(k)+(t_end-1)*dt, 0 );

            [ ~ , X_eta ] = Wingbeat_Legendre_matrix( n_pol_eta, down_up_ratio(k), 100, time(k), time(k)+(t_end-1)*dt, 0 );

            [ t_temp , X_phi ] = Wingbeat_Legendre_matrix( n_pol_phi, down_up_ratio(k), 100, time(k), time(k)+(t_end-1)*dt, 0 );

            time2(((k-1)*100+1):(k*100)) = t_temp;

            [ ~ , X_dot_theta ] = Wingbeat_Legendre_matrix( n_pol_theta, down_up_ratio(k), 100, time(k), time(k)+(t_end-1)*dt, 1 );

            [ ~ , X_dot_eta ] = Wingbeat_Legendre_matrix( n_pol_eta, down_up_ratio(k), 100, time(k), time(k)+(t_end-1)*dt, 1 );

            [ ~ , X_dot_phi ] = Wingbeat_Legendre_matrix( n_pol_phi, down_up_ratio(k), 100, time(k), time(k)+(t_end-1)*dt, 1 );

            dev_theta_L(((k-1)*100+1):(k*100)) = X_theta*[a_dev.theta_L1(:,k); a_dev.theta_L2(:,k)];

            dev_theta_R(((k-1)*100+1):(k*100)) = X_theta*[a_dev.theta_R1(:,k); a_dev.theta_R2(:,k)];

            dev_eta_L(((k-1)*100+1):(k*100)) = X_eta*[a_dev.eta_L1(:,k); a_dev.eta_L2(:,k)];

            dev_eta_R(((k-1)*100+1):(k*100)) = X_eta*[a_dev.eta_R1(:,k); a_dev.eta_R2(:,k)];

            dev_phi_L(((k-1)*100+1):(k*100)) = X_phi*[a_dev.phi_L1(:,k); a_dev.phi_L2(:,k)];

            dev_phi_R(((k-1)*100+1):(k*100)) = X_phi*[a_dev.phi_R1(:,k); a_dev.phi_R2(:,k)];

            dev_theta_dot_L(((k-1)*100+1):(k*100)) = X_dot_theta*[a_dev.theta_L1(:,k); a_dev.theta_L2(:,k)];

            dev_theta_dot_R(((k-1)*100+1):(k*100)) = X_dot_theta*[a_dev.theta_R1(:,k); a_dev.theta_R2(:,k)];

            dev_eta_dot_L(((k-1)*100+1):(k*100)) = X_dot_eta*[a_dev.eta_L1(:,k); a_dev.eta_L2(:,k)];

            dev_eta_dot_R(((k-1)*100+1):(k*100)) = X_dot_eta*[a_dev.eta_R1(:,k); a_dev.eta_R2(:,k)];

            dev_phi_dot_L(((k-1)*100+1):(k*100)) = X_dot_phi*[a_dev.phi_L1(:,k); a_dev.phi_L2(:,k)];

            dev_phi_dot_R(((k-1)*100+1):(k*100)) = X_dot_phi*[a_dev.phi_R1(:,k); a_dev.phi_R2(:,k)];    
            
        end

    end
    
    dev_data = {};
    
    dev_data.theta_L = dev_theta_L;
    
    dev_data.theta_R = dev_theta_R;
    
    dev_data.eta_L = dev_eta_L;
    
    dev_data.eta_R = dev_eta_R;
    
    dev_data.phi_L = dev_phi_L;
    
    dev_data.phi_R = dev_phi_R;
    
    
    ABCD = ABCD_points( a_fit, a_avg, down_up_ratio, down_up_avg, n_pol_theta, n_pol_eta, n_pol_phi, nr_wb );
    
    if i == 1
    
        A_theta_L = ABCD.dev_L_theta(:,1);
        A_eta_L = ABCD.dev_L_eta(:,1);
        A_phi_L = ABCD.dev_L_phi(:,1);

        B_theta_L = ABCD.dev_L_theta(:,2);
        B_eta_L = ABCD.dev_L_eta(:,2);
        B_phi_L = ABCD.dev_L_phi(:,2);

        C_theta_L = ABCD.dev_L_theta(:,3);
        C_eta_L = ABCD.dev_L_eta(:,3);
        C_phi_L = ABCD.dev_L_phi(:,3);

        D_theta_L = ABCD.dev_L_theta(:,4);
        D_eta_L = ABCD.dev_L_eta(:,4);
        D_phi_L = ABCD.dev_L_phi(:,4);
    
        A_theta_R = ABCD.dev_R_theta(:,1);
        A_eta_R = ABCD.dev_R_eta(:,1);
        A_phi_R = ABCD.dev_R_phi(:,1);

        B_theta_R= ABCD.dev_R_theta(:,2);
        B_eta_R = ABCD.dev_R_eta(:,2);
        B_phi_R = ABCD.dev_R_phi(:,2);

        C_theta_R = ABCD.dev_R_theta(:,3);
        C_eta_R = ABCD.dev_R_eta(:,3);
        C_phi_R = ABCD.dev_R_phi(:,3);

        D_theta_R = ABCD.dev_R_theta(:,4);
        D_eta_R = ABCD.dev_R_eta(:,4);
        D_phi_R = ABCD.dev_R_phi(:,4);
        
    else
        
        A_theta_L = [A_theta_L; ABCD.dev_L_theta(:,1)];
        A_eta_L = [A_eta_L; ABCD.dev_L_eta(:,1)];
        A_phi_L = [A_phi_L; ABCD.dev_L_phi(:,1)];

        B_theta_L = [B_theta_L; ABCD.dev_L_theta(:,2)];
        B_eta_L = [B_eta_L; ABCD.dev_L_eta(:,2)];
        B_phi_L = [B_phi_L; ABCD.dev_L_phi(:,2)];

        C_theta_L = [C_theta_L; ABCD.dev_L_theta(:,3)];
        C_eta_L = [C_eta_L; ABCD.dev_L_eta(:,3)];
        C_phi_L = [C_phi_L; ABCD.dev_L_phi(:,3)];

        D_theta_L = [D_theta_L; ABCD.dev_L_theta(:,4)];
        D_eta_L = [D_eta_L; ABCD.dev_L_eta(:,4)];
        D_phi_L = [D_phi_L; ABCD.dev_L_phi(:,4)];
    
        A_theta_R = [A_theta_R; ABCD.dev_R_theta(:,1)];
        A_eta_R = [A_eta_R; ABCD.dev_R_eta(:,1)];
        A_phi_R = [A_phi_R; ABCD.dev_R_phi(:,1)];

        B_theta_R= [B_theta_R; ABCD.dev_R_theta(:,2)];
        B_eta_R = [B_eta_R; ABCD.dev_R_eta(:,2)];
        B_phi_R = [B_phi_R; ABCD.dev_R_phi(:,2)];

        C_theta_R = [C_theta_R; ABCD.dev_R_theta(:,3)];
        C_eta_R = [C_eta_R; ABCD.dev_R_eta(:,3)];
        C_phi_R = [C_phi_R; ABCD.dev_R_phi(:,3)];

        D_theta_R = [D_theta_R; ABCD.dev_R_theta(:,4)];
        D_eta_R = [D_eta_R; ABCD.dev_R_eta(:,4)];
        D_phi_R = [D_phi_R; ABCD.dev_R_phi(:,4)];      
    
    end

    
    c_var = 1.5; % threshold maneuvering flight ( c_var times variance steady flight ).
    
    maneuver_wb = select_maneuver( dev_data,trigger_wb, nr_wb, c_var );


    maneuver_wb.id
    
    
    nr_of_maneuvers = length(maneuver_wb.id);
    
    temp_maneuver_matrix = zeros(nr_of_maneuvers,25);
    
    for j = 1:nr_of_maneuvers
        
        temp_maneuver_matrix(j,:) = [ maneuver_wb.id(j) i down_up_ratio(j)  stroke_var.Omega_strk(:,maneuver_wb.id(j))'     stroke_var.Omega_dot_strk(:,maneuver_wb.id(j))' ...
                                      stroke_var.Vn(maneuver_wb.id(j))      stroke_var.Vt(maneuver_wb.id(j))        stroke_var.An(maneuver_wb.id(j))        stroke_var.At(maneuver_wb.id(j)) ...
                                      ABCD.dev_L_theta(maneuver_wb.id(j),1) ABCD.dev_L_eta(maneuver_wb.id(j),1)     ABCD.dev_L_phi(maneuver_wb.id(j),1)     ABCD.dev_L_theta(maneuver_wb.id(j),2) ...
                                      ABCD.dev_L_eta(maneuver_wb.id(j),2)   ABCD.dev_L_phi(maneuver_wb.id(j),2)     ABCD.dev_L_theta(maneuver_wb.id(j),3)   ABCD.dev_L_eta(maneuver_wb.id(j),3) ...
                                      ABCD.dev_L_phi(maneuver_wb.id(j),3)   ABCD.dev_L_theta(maneuver_wb.id(j),4)   ABCD.dev_L_eta(maneuver_wb.id(j),4)     ABCD.dev_L_phi(maneuver_wb.id(j),4) ];
        
    end
    
    
    
  
    
    
    if i == 1
        
        Maneuver_Matrix = temp_maneuver_matrix;
        
        a_sym_maneuver.theta_L = [a_sym.theta_L1; a_sym.theta_L2];

        a_sym_maneuver.theta_R = [a_sym.theta_R1; a_sym.theta_R2];

        a_sym_maneuver.eta_L = [a_sym.eta_L1; a_sym.eta_L2];

        a_sym_maneuver.eta_R = [a_sym.eta_R1; a_sym.eta_R2];

        a_sym_maneuver.phi_L = [a_sym.phi_L1; a_sym.phi_L2];

        a_sym_maneuver.phi_R = [a_sym.phi_R1; a_sym.phi_R2];
       
    
        a_dev_maneuver.theta_L = [a_dev.theta_L1; a_dev.theta_L2];

        a_dev_maneuver.theta_R = [a_dev.theta_R1; a_dev.theta_R2];

        a_dev_maneuver.eta_L = [a_dev.eta_L1; a_dev.eta_L2];

        a_dev_maneuver.eta_R = [a_dev.eta_R1; a_dev.eta_R2];

        a_dev_maneuver.phi_L = [a_dev.phi_L1; a_dev.phi_L2];

        a_dev_maneuver.phi_R = [a_dev.phi_R1; a_dev.phi_R2];
        
    else
        
        Maneuver_Matrix = [Maneuver_Matrix; temp_maneuver_matrix];
        
        a_sym_maneuver.theta_L = [a_sym_maneuver.theta_L [a_sym.theta_L1; a_sym.theta_L2]];

        a_sym_maneuver.theta_R = [a_sym_maneuver.theta_R [a_sym.theta_R1; a_sym.theta_R2]];

        a_sym_maneuver.eta_L = [a_sym_maneuver.eta_L [a_sym.eta_L1; a_sym.eta_L2]];

        a_sym_maneuver.eta_R = [a_sym_maneuver.eta_R [a_sym.eta_R1; a_sym.eta_R2]];

        a_sym_maneuver.phi_L = [a_sym_maneuver.phi_L [a_sym.phi_L1; a_sym.phi_L2]];

        a_sym_maneuver.phi_R = [a_sym_maneuver.phi_R [a_sym.phi_R1; a_sym.phi_R2]];
       
    
        a_dev_maneuver.theta_L = [a_dev_maneuver.theta_L [a_dev.theta_L1; a_dev.theta_L2]];

        a_dev_maneuver.theta_R = [a_dev_maneuver.theta_R [a_dev.theta_R1; a_dev.theta_R2]];

        a_dev_maneuver.eta_L = [a_dev_maneuver.eta_L [a_dev.eta_L1; a_dev.eta_L2]];

        a_dev_maneuver.eta_R = [a_dev_maneuver.eta_R [a_dev.eta_R1; a_dev.eta_R2]];

        a_dev_maneuver.phi_L = [a_dev_maneuver.phi_L [a_dev.phi_L1; a_dev.phi_L2]];

        a_dev_maneuver.phi_R = [a_dev_maneuver.phi_R [a_dev.phi_R1; a_dev.phi_R2]];
        
    end
    
    
    
    
%     figure()
%     hold on
%     subplot(3,1,1); plot(time2,radtodeg(dev_theta_L),'r',time2,radtodeg(dev_theta_R),'g');
%     subplot(3,1,2); plot(time2,radtodeg(dev_eta_L),'r',time2,radtodeg(dev_eta_R),'g');
%     subplot(3,1,3); plot(time2,radtodeg(dev_phi_L),'r',time2,radtodeg(dev_phi_R),'g');
%     hold off
% 
%     figure()
%     hold on
%     subplot(2,1,1); plot(time2,radtodeg(dev_eta_L),'r',time2,radtodeg(dev_eta_R),'g');
%     subplot(2,1,2); plot(time2,radtodeg(sqrt(dev_theta_L.^2+dev_phi_L.^2)),'r',time2,radtodeg(sqrt(dev_theta_R.^2+dev_phi_R.^2)),'g')
% 
%     figure()
%     hold on
%     subplot(3,1,1); plot(stroke_var.Omega_strk(1,:))
%     subplot(3,1,2); plot(stroke_var.Omega_strk(2,:))
%     subplot(3,1,3); plot(stroke_var.Omega_strk(3,:))
%     hold off
%     
%     pause

    

    end
    
    end

    
    a_sym_maneuver
    
    a_dev_maneuver


    Maneuver_Matrix


    
%     figure()
%     plot(A_phi_L,A_theta_L,'o')
%     
%     figure()
%     plot(B_phi_L,B_theta_L,'o')
%     
%     figure()
%     plot(C_phi_L,C_theta_L,'o')
%     
%     figure()
%     plot(D_phi_L,D_theta_L,'o')
%     
%     figure()
%     plot(A_phi_R,A_theta_R,'o')
%     
%     figure()
%     plot(B_phi_R,B_theta_R,'o')
%     
%     figure()
%     plot(C_phi_R,C_theta_R,'o')
%     
%     figure()
%     plot(D_phi_R,D_theta_R,'o')
    
    figure()
    plot3(Maneuver_Matrix(:,14),Maneuver_Matrix(:,16),Maneuver_Matrix(:,7),'o')
    title('point A over strokeplane roll acceleration')
    xlabel('theta')
    ylabel('phi')
    zlabel('omega x dot')
    
    figure()
    plot3(Maneuver_Matrix(:,17),Maneuver_Matrix(:,19),Maneuver_Matrix(:,7),'o')
    title('point B over strokeplane roll acceleration')
    xlabel('theta')
    ylabel('phi')
    zlabel('omega x dot')
    
    figure()
    plot3(Maneuver_Matrix(:,20),Maneuver_Matrix(:,22),Maneuver_Matrix(:,7),'o')
    title('point C over strokeplane roll acceleration')
    xlabel('theta')
    ylabel('phi')
    zlabel('omega x dot')

    figure()
    plot3(Maneuver_Matrix(:,23),Maneuver_Matrix(:,25),Maneuver_Matrix(:,7),'o')
    title('point D over strokeplane roll acceleration')
    xlabel('theta')
    ylabel('phi')
    zlabel('omega x dot')
    
    figure()
    plot3(Maneuver_Matrix(:,20)-Maneuver_Matrix(:,23),Maneuver_Matrix(:,22)-Maneuver_Matrix(:,25),Maneuver_Matrix(:,7),'o')
    title('point C-D over strokeplane roll acceleration')
    xlabel('theta')
    ylabel('phi')
    zlabel('omega x dot')

    [ ABCD_fit ] = ABCD_fitting( Maneuver_Matrix, dt )


