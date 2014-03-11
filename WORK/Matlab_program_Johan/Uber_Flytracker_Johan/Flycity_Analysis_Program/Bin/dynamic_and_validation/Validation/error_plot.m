function error_plot( x_strk_sim_steady, x_strk_mov_steady, x_strk_sim_man, x_strk_mov_man, man_id, i)

    % Compute the error vector and the error angle
    
    N_steady = size(x_strk_sim_steady,2);
    
    for j = 1:N_steady
   
%     err_vect_steady(j)   = norm(x_strk_sim_steady(:,j)-x_strk_mov_steady(:,j));
%     err_angle_steady(j)  = norm(radtodeg(asin(norm(cross(x_strk_sim_steady(:,j),x_strk_mov_steady(:,j)))/(norm(x_strk_sim_steady(:,j))*norm(x_strk_mov_steady(:,j))))));

    err_vect_steady(j)   = norm(x_strk_sim_steady(i,j)-x_strk_mov_steady(i,j));
    
    if i == 1
        
        err_angle_steady(j)  = radtodeg(norm(atan(x_strk_sim_steady(3,j)/x_strk_sim_steady(2,j))-atan(x_strk_mov_steady(3,j)/x_strk_mov_steady(2,j))));
        
    elseif i == 2
        
        err_angle_steady(j)  = radtodeg(norm(atan(x_strk_sim_steady(3,j)/x_strk_sim_steady(1,j))-atan(x_strk_mov_steady(3,j)/x_strk_mov_steady(1,j))));
        
    elseif i == 3
        
        err_angle_steady(j)  = radtodeg(norm(atan(x_strk_sim_steady(2,j)/x_strk_sim_steady(1,j))-atan(x_strk_mov_steady(2,j)/x_strk_mov_steady(1,j))));
        
    end
    
    end
    
    N_man = size(x_strk_sim_man,2);
    
    for j = 1:N_man
    
%     err_vect_man(j)      = norm(x_strk_sim_man(:,j)-x_strk_mov_man(:,j));
%     err_angle_man(j)     = norm(radtodeg(asin(norm(cross(x_strk_sim_man(:,j),x_strk_mov_man(:,j)))/(norm(x_strk_sim_man(:,j))*norm(x_strk_mov_man(:,j))))));

    err_vect_man(j)      = norm(x_strk_sim_man(i,j)-x_strk_mov_man(i,j));
    
    if i == 1
        
        err_angle_man(j)  = radtodeg(norm(atan(x_strk_sim_man(3,j)/x_strk_sim_man(2,j))-atan(x_strk_mov_man(3,j)/x_strk_mov_man(2,j))));
        
    elseif i == 2
        
        err_angle_man(j)  = radtodeg(norm(atan(x_strk_sim_man(3,j)/x_strk_sim_man(1,j))-atan(x_strk_mov_man(3,j)/x_strk_mov_man(1,j))));
        
    elseif i == 3
        
        err_angle_man(j)  = radtodeg(norm(atan(x_strk_sim_man(2,j)/x_strk_sim_man(1,j))-atan(x_strk_mov_man(2,j)/x_strk_mov_man(1,j))));
        
    end

    end
    
    % Compute the average and the std:
    
    n_bins = 9;
    
    x_mov_all           = [ x_strk_mov_steady(i,:) x_strk_mov_man(i,:) ];
    
    err_vect_all        = [ err_vect_steady err_vect_man ];
    err_angle_all       = [ err_angle_steady err_angle_man ];

    [x_mov_all_sort, sort_id_x] = sort(x_mov_all);
    
    err_vect_all_sort   = err_vect_all(sort_id_x);
    err_angle_all_sort  = err_angle_all(sort_id_x);

    N = length(err_angle_all_sort);
    
    dN = (N-1)/n_bins;
    
    bin_loc             = zeros(n_bins,1);
    
    bin_mean_vect       = zeros(n_bins,1);
    bin_std_vect        = zeros(n_bins,1);
    
    bin_mean_angle      = zeros(n_bins,1);
    bin_std_angle       = zeros(n_bins,1);
    
    N_loc = round(1:dN:N);
        
    for k = 1:n_bins
        
        t_range = N_loc(k):N_loc(k+1);
        
        bin_loc(k)          = mean(x_mov_all_sort(t_range));
        
        [ w_avg_vect, ~, w_std_vect ] = weighted_avg(err_vect_all_sort(t_range));
        
        bin_mean_vect(k)    = w_avg_vect;
        bin_std_vect(k)     = w_std_vect;
        
        [ w_avg_angle, ~, w_std_angle ] = weighted_avg(err_angle_all_sort(t_range));
        
        bin_mean_angle(k)   = w_avg_angle;
        bin_std_angle(k)    = w_std_angle;
        
    end
    
    % Plot:
    
    figure()
    hold on
    scatter(x_strk_mov_steady(i,:),err_vect_steady,'.','r')
    scatter(x_strk_mov_man(i,:),err_vect_man,'.','b')
    plot(bin_loc,bin_mean_vect,'k')
    plot(bin_loc,bin_mean_vect+bin_std_vect,'k')
    plot(bin_loc,bin_mean_vect-bin_std_vect,'k')
    hold off
    
    figure()
    hold on
    scatter(x_strk_mov_steady(i,:),err_angle_steady,'.','r')
    scatter(x_strk_mov_man(i,:),err_angle_man,'.','b')
    plot(bin_loc,bin_mean_angle,'k')
    plot(bin_loc,bin_mean_angle+bin_std_angle,'k')
    plot(bin_loc,bin_mean_angle-bin_std_angle,'k')
    hold off
    
    
end