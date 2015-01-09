function Legendre_Bernstein_fitting( settings, pathDB)



    % This program fits polynomial functions to the wingbeat kinematics
    % according to the CST geometry parameterization. This CST
    % representation of class / shape transformations combines the
    % strengths of two classes of polynomials. The Legendre polynomials
    % have their strength in the general representation of the wingbeat as
    % they can represent harmonic and near-harmonic functions rather well.
    % A disadvantage of the Legendre polynomials is the introduction of
    % higher order oscillations in the approximation when the number of
    % polynomials is large. This is solved in the CST representation by
    % taking a low order Legendre representation for the class functions
    % and a local refinement of this class function by Bernstein
    % polynomials in the shape function. In this way the polynomial fit
    % doesn't introduce unwanted oscillations whilst the parameterization
    % of the polynomial fit remains small.
    
    
    % Loop:
    
    for i =1:size(pathDB.x,2)
        
    
    % Initialization:
    
    seq_nr = i

    % Nr Legendre polynomials:
   
    n_Leg_pol_theta = 8; % Order of used Legendre polynomials
    n_Leg_pol_eta = 8; % Order of used Legendre polynomials
    n_Leg_pol_phi = 8; % Order of used Legendre polynomials
    
    
    n_pol_Leg = {};
    
    n_pol_Leg.theta = n_Leg_pol_theta;
    
    n_pol_Leg.eta = n_Leg_pol_eta;
    
    n_pol_Leg.phi = n_Leg_pol_phi;
    
    % Nr Bernstein polynomials:
    
    n_Bern_pol_theta = 8; % Order of used Bernstein polynomials
    n_Bern_pol_eta = 8; % Order of used Bernstein polynomials
    n_Bern_pol_phi = 8; % Order of used Bernstein polynomials
    
    
    n_pol_Bern = {};
    
    n_pol_Bern.theta = n_Bern_pol_theta;
    
    n_pol_Bern.eta = n_Bern_pol_eta;
    
    n_pol_Bern.phi = n_Bern_pol_phi;    
   
    
    % Wing kinematic angles:
    
    start = find(isnan(pathDB.x(:,seq_nr))==0, 1 );
    stop = find(isnan(pathDB.x(:,seq_nr))==0, 1, 'last' );
    
    dt = pathDB.t(2)-pathDB.t(1);

    eta_L = pathDB.eta_L(start:stop,seq_nr);

    theta_L = pathDB.theta_L(start:stop,seq_nr);

    phi_L = pathDB.phi_L(start:stop,seq_nr);

    eta_R = pathDB.eta_R(start:stop,seq_nr);

    theta_R = pathDB.theta_R(start:stop,seq_nr);

    phi_R = pathDB.phi_R(start:stop,seq_nr);
    
    
    % Nr of wingbeats:

    nr_wb = find(isnan(pathDB.wingbeat_time(:,1,seq_nr))==0, 1, 'last' );
    
    
    % Wingbeat location:

    wb_loc = zeros(nr_wb,2);

    for k = 1:nr_wb

        if k == nr_wb

            wb_end = find(isnan(pathDB.wingbeat_time(nr_wb,:,seq_nr))==0, 1, 'last' );

            wb_loc(k,:) = [ pathDB.wingbeat_time(nr_wb,1,seq_nr) pathDB.wingbeat_time(nr_wb,wb_end,seq_nr) ];

        else

            wb_loc(k,:) = [ pathDB.wingbeat_time(k,1,seq_nr) pathDB.wingbeat_time(k+1,1,seq_nr) ];

        end

    end
    
    
    % Legendre fit coefficients:
    
    a_fit_Leg = {};
    
    
    % Bernstein fit coefficients:
    
    a_fit_Bern = {};
    
    % Step 1: Class function fit.
    
    a_fit_Leg.eta_L = Piecewise_polynomial_fit(eta_L,n_pol_Leg.eta,wb_loc);
    
    a_fit_Leg.theta_L = Piecewise_polynomial_fit(theta_L,n_pol_Leg.theta,wb_loc);
    
    a_fit_Leg.phi_L = Piecewise_polynomial_fit(phi_L,n_pol_Leg.phi,wb_loc);
    
    a_fit_Leg.eta_R = Piecewise_polynomial_fit(eta_R,n_pol_Leg.eta,wb_loc);
    
    a_fit_Leg.theta_R = Piecewise_polynomial_fit(theta_R,n_pol_Leg.theta,wb_loc);
    
    a_fit_Leg.phi_R = Piecewise_polynomial_fit(phi_R,n_pol_Leg.phi,wb_loc);
    
    
    % Step 2: Shape function fit.
    
    
    
    % Step 3: Determine standard wingbeat and maneuver wingbeat.
    
    
    
    % End loop.
    
    end
    
    % Step 4: Output.


end

