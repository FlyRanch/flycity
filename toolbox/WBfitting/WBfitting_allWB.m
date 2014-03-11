function [a_fit, a_avg] = WBfitting_allWB(n_pol_stroke,n_pol_pitch,n_pol_dev,...
    t_loc,stroke_loc,pitch_loc,dev_loc)


a_fit = {};

n_start = find(t_loc == 0);
n_stop = find(t_loc == 1);
N_wb = length(n_stop);

    for j = 1:N_wb
        
        [ a_stroke1, a_stroke2, ~ ] = local_fit(stroke_loc(n_start(j):n_stop(j)),n_pol_stroke);
        
        [ a_pitch1, a_pitch2, ~ ] = local_fit(pitch_loc(n_start(j):n_stop(j)),n_pol_pitch);
        
        [ a_dev1, a_dev2, ratio_1_2 ] = local_fit(dev_loc(n_start(j):n_stop(j)),n_pol_dev);
        
        a_fit_stroke1(j,:) = a_stroke1;
        
        a_fit_pitch1(j,:) = a_pitch1;
        
        a_fit_dev1(j,:) = a_dev1;
        
        a_fit_stroke2(j,:) = a_stroke2;
        
        a_fit_pitch2(j,:) = a_pitch2;
        
        a_fit_dev2(j,:) = a_dev2;
        
        a_fit_ratio_1_2(j) = ratio_1_2;
        
        a_fit_wb1_wb2(j,1) = ceil((ratio_1_2*(n_stop(j)-n_start(j)+1))/2);
        
        a_fit_wb1_wb2(j,2) = floor(((2-ratio_1_2)*(n_stop(j)-n_start(j)+1))/2);
        
    end    
    
            a_fit.stroke1 = a_fit_stroke1;
            a_fit.stroke2 = a_fit_stroke2;
            a_fit.pitch1 = a_fit_pitch1;
            a_fit.pitch2 = a_fit_pitch2;
            a_fit.dev1 = a_fit_dev1;
            a_fit.dev2 = a_fit_dev2;

            a_fit.ratio_1_2 = a_fit_ratio_1_2;
            
            a_fit.wb1_wb2 = a_fit_wb1_wb2;
            
    
    
    % Adjust wb_loc_mean and stroke etc. then it should work.
        
    wb_loc_mean12 = zeros(2*length(a_fit.wb1_wb2(:,1)),2);
    
    for k = 1:(length(a_fit.wb1_wb2(:,1)))
               
        wb_loc_mean12(2*k-1,:) = [1 a_fit.wb1_wb2(k,1)];
        
        wb_loc_mean12(2*k,:) = [a_fit.wb1_wb2(k,1)+1 (a_fit.wb1_wb2(k,1)+a_fit.wb1_wb2(k,2))];

    end
        
        [a_avg_stroke1, a_avg_stroke2] = average_fit(a_fit.stroke1,a_fit.stroke2,n_pol_stroke);
        
        [a_avg_pitch1, a_avg_pitch2] = average_fit(a_fit.pitch1,a_fit.pitch2,n_pol_pitch);
        
        [a_avg_dev1, a_avg_dev2] = average_fit(a_fit.dev1,a_fit.dev2,n_pol_dev);
        
        % Store average coefficients in a_avg
        
        a_avg = {};

        a_avg.stroke1 = a_avg_stroke1;
        a_avg.stroke2 = a_avg_stroke2;
        a_avg.pitch1 = a_avg_pitch1;
        a_avg.pitch2 = a_avg_pitch2;
        a_avg.dev1 = a_avg_dev1;
        a_avg.dev2 = a_avg_dev2;
        
        a_avg.wb1_wb2 = [mean(a_fit.wb1_wb2(:,1)) mean(a_fit.wb1_wb2(:,2))];

    
    
    
    % Plot results:
    
    PN_stroke = Legendre_polynomial(n_pol_stroke,2,-1:0.01:1);
    PN_pitch = Legendre_polynomial(n_pol_pitch,2,-1:0.01:1);
    PN_dev = Legendre_polynomial(n_pol_dev,2,-1:0.01:1);

    nr_wb = length(a_fit.wb1_wb2(:,1))
    
%     figure()
    subplot(3,5,1)
    hold on
%     for k = 1:nr_wb
%         m1 =a_fit.wb1_wb2(k,1);
%         m2 = a_fit.wb1_wb2(k,2);
%         t1 = -1:((m1/m2)/200):(-1+(m1/m2));
%         t2 = (-1+(m1/m2)):((2-(m1/m2))/200):1;
%         t1_mod=(t1+1)/2;
%         t2_mod=(t2+1)/2;
%         plot(t1_mod,(PN_stroke(:,:,1)'*a_fit.stroke1(k,:)'),'Color',[0.5 0.5 0.5])
%         plot(t2_mod,(PN_stroke(:,:,1)'*a_fit.stroke2(k,:)'),'Color',[0.5 0.5 0.5])
%     end
    m1 =a_avg.wb1_wb2(1);
    m2 = a_avg.wb1_wb2(2);
    t1 = -1:((m1/m2)/200):(-1+(m1/m2));
    t2 = (-1+(m1/m2)):((2-(m1/m2))/200):1;
        t1_mod=(t1+1)/2;
        t2_mod=(t2+1)/2;
    plot(t1_mod,(PN_stroke(:,:,1)'*a_avg.stroke1),'r')
    plot(t2_mod,(PN_stroke(:,:,1)'*a_avg.stroke2),'r')    
    hold off
    
%     figure()
    subplot(3,5,2)
    hold on
%     for k = 1:nr_wb
%         m1 =a_fit.wb1_wb2(k,1);
%         m2 = a_fit.wb1_wb2(k,2);
%         t1 = -1:((m1/m2)/200):(-1+(m1/m2));
%         t2 = (-1+(m1/m2)):((2-(m1/m2))/200):1;
%         t1_mod=(t1+1)/2;
%         t2_mod=(t2+1)/2;
%         plot(t1_mod,(PN_pitch(:,:,1)'*a_fit.pitch1(k,:)'),'Color',[0.5 0.5 0.5])
%         plot(t2_mod,(PN_pitch(:,:,1)'*a_fit.pitch2(k,:)'),'Color',[0.5 0.5 0.5])
%     end
    m1 =a_avg.wb1_wb2(1);
    m2 = a_avg.wb1_wb2(2);
    t1 = -1:((m1/m2)/200):(-1+(m1/m2));
    t2 = (-1+(m1/m2)):((2-(m1/m2))/200):1;
        t1_mod=(t1+1)/2;
        t2_mod=(t2+1)/2;
    plot(t1_mod,(PN_pitch(:,:,1)'*a_avg.pitch1),'r')
    plot(t2_mod,(PN_pitch(:,:,1)'*a_avg.pitch2),'r')   
    hold off
    
%     figure()
    subplot(3,5,3)
    hold on
%     for k = 1:nr_wb
%         m1 =a_fit.wb1_wb2(k,1);
%         m2 = a_fit.wb1_wb2(k,2);
%         t1 = -1:((m1/m2)/200):(-1+(m1/m2));
%         t2 = (-1+(m1/m2)):((2-(m1/m2))/200):1;
%         t1_mod=(t1+1)/2;
%         t2_mod=(t2+1)/2;
%         plot(-1:0.005:0,(PN_dev(:,:,1)'*a_fit.dev1(k,:)'),'Color',[0.5 0.5 0.5])
%         plot(0:0.005:1,(PN_dev(:,:,1)'*a_fit.dev2(k,:)'),'Color',[0.5 0.5 0.5])
%     end
    m1 =a_avg.wb1_wb2(1);
    m2 = a_avg.wb1_wb2(2);
    t1 = -1:((m1/m2)/200):(-1+(m1/m2));
    t2 = (-1+(m1/m2)):((2-(m1/m2))/200):1;
        t1_mod=(t1+1)/2;
        t2_mod=(t2+1)/2;
    plot(t1_mod,(PN_dev(:,:,1)'*a_avg.dev1),'r')
    plot(t2_mod,(PN_dev(:,:,1)'*a_avg.dev2),'r')   
    hold off
    

