function pathDB8( settings, pathDB )


    savefile = 'pathDB8.mat';

    n_pol_theta = 12; % Order of used polynomials
    n_pol_eta = 14; % Order of used polynomials
    n_pol_phi = 10; % Order of used polynomials
    
    w_tresh = 6*pi;
    
    a_tresh = 4000;
    
        
    dt = pathDB.t(2)-pathDB.t(1);
    
    nr_of_seq = size(pathDB.x,2);   

    
    a_fit_tot = [];
    a_avg_tot = [];
    a_sym_tot = [];
    a_dev_tot = [];
    f_avg_tot = [];
    down_up_avg_tot = [];
    trigger_wb_tot = [];
    down_up_ratio_tot = [];
    stroke_var_tot = [];
    maneuver_tot = [];
    
    for i = 1:nr_of_seq
        
        if i == 1
            
            a_fit_tot = char(['a_fit_' int2str(i)]);
            a_avg_tot = char(['a_avg_' int2str(i)]);
            a_sym_tot = char(['a_sym_' int2str(i)]);
            a_dev_tot = char(['a_dev_' int2str(i)]);
            f_avg_tot = char(['f_avg_' int2str(i)]);
            down_up_avg_tot = char(['down_up_avg_' int2str(i)]);
            trigger_wb_tot = char(['trigger_wb_' int2str(i)]);
            down_up_ratio_tot = char(['down_up_ratio_' int2str(i)]);
            stroke_var_tot = char(['stroke_var_' int2str(i)]);
            maneuver_tot = char(['maneuver_' int2str(i)]);
            
        else
            
            a_fit_tot = [a_fit_tot; char(['a_fit_' int2str(i)])];
            a_avg_tot = [a_avg_tot; char(['a_avg_' int2str(i)])];
            a_sym_tot = [a_sym_tot; char(['a_sym_' int2str(i)])];
            a_dev_tot = [a_dev_tot; char(['a_dev_' int2str(i)])];
            f_avg_tot = [f_avg_tot; char(['f_avg_' int2str(i)])];
            down_up_avg_tot = [down_up_avg_tot; char(['down_up_avg_' int2str(i)])];
            trigger_wb_tot = [trigger_wb_tot; char(['trigger_wb_' int2str(i)])];
            down_up_ratio_tot = [down_up_ratio_tot; char(['down_up_ratio_' int2str(i)])];
            stroke_var_tot = [stroke_var_tot; char(['stroke_var_' int2str(i)])];
            maneuver_tot = [maneuver_tot; char(['maneuver_' int2str(i)])];            
            
        end
        
    end
    
    
    for i = 1:nr_of_seq
        
        seq_nr = i

%         start = find(isnan(pathDB.x(:,seq_nr))==0, 1 );
%         stop = find(isnan(pathDB.x(:,seq_nr))==0, 1, 'last' );

        nr_wb = find(isnan(pathDB.wingbeat_time(:,1,seq_nr))==0, 1, 'last' );
        
        [a_fit,a_avg,f_avg,down_up_avg,trigger_wb,down_up_ratio] = Standard_wingbeat( settings, pathDB, i, n_pol_theta, n_pol_eta, n_pol_phi);
        
        [a_sym, a_dev] = Maneuver_wingbeat(a_fit,a_avg,n_pol_theta,n_pol_eta,n_pol_phi,down_up_ratio,down_up_avg);
        
        [ stroke_var ] = Strokeplane_dynamics( settings, pathDB, i );
        
        [ maneuver ] = select_maneuvers2(stroke_var,w_tresh,a_tresh);
                    
        Maneuver_matrix = [ Maneuver_matrix; ...
                           [maneuver.roll_turns' maneuver.pitch_turns' maneuver.yaw_turns' maneuver.a_n' maneuver.a_t' [1:nr_wb]' ones(nr_wb,1)*i] ];
                       
        a_fit_tot(i) = a_fit;
        a_avg_tot(i) = a_avg;
        a_sym_tot(i) = a_sym;
        a_dev_tot(i) = a_dev;
        f_avg_tot(i) = f_avg;
        down_up_avg_tot(i) = down_up_avg;
        trigger_wb_tot(i) = trigger_wb;
        down_up_ratio_tot(i) = down_up_ratio;
        stroke_var_tot(i) = stroke_var;
        maneuver_tot(i) = maneuver;
                            
    end



    save(savefile,'a_fit_tot','a_avg_tot','a_sym_tot','a_dev_tot', 'f_avg_tot', 'down_up_avg_tot', 'trigger_wb_tot', 'down_up_ratio_tot', 'stroke_var_tot', 'maneuver_tot' )

end

