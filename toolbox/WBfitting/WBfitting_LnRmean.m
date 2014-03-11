function [a_fit, a_avg] = WBfitting_LnRmean(n_pol_stroke,n_pol_pitch,n_pol_dev,...
    t_loc,stroke_L_loc,pitch_L_loc,dev_L_loc,...
    stroke_R_loc,pitch_R_loc,dev_R_loc)


a_fit = {};

n_start = find(t_loc == 0);
n_stop = find(t_loc == 1);
N_wb = length(n_stop);

    for j = 1:N_wb
        
        [ a_stroke_L1, a_stroke_L2, ~ ] = local_fit(stroke_L_loc(n_start(j):n_stop(j)),n_pol_stroke);
        
        [ a_pitch_L1, a_pitch_L2, ~ ] = local_fit(pitch_L_loc(n_start(j):n_stop(j)),n_pol_pitch);
        
        [ a_dev_L1, a_dev_L2, ~ ] = local_fit(dev_L_loc(n_start(j):n_stop(j)),n_pol_dev);

        [ a_stroke_R1, a_stroke_R2, ~ ] = local_fit(stroke_R_loc(n_start(j):n_stop(j)),n_pol_stroke);
        
        [ a_pitch_R1, a_pitch_R2, ~ ] = local_fit(pitch_R_loc(n_start(j):n_stop(j)),n_pol_pitch);
        
        [ a_dev_R1, a_dev_R2, ratio_1_2 ] = local_fit(dev_R_loc(n_start(j):n_stop(j)),n_pol_dev);
        
        a_fit_stroke_L1(j,:) = a_stroke_L1;
        
        a_fit_pitch_L1(j,:) = a_pitch_L1;
        
        a_fit_dev_L1(j,:) = a_dev_L1;
        
        a_fit_stroke_L2(j,:) = a_stroke_L2;
        
        a_fit_pitch_L2(j,:) = a_pitch_L2;
        
        a_fit_dev_L2(j,:) = a_dev_L2;

        a_fit_stroke_R1(j,:) = a_stroke_R1;
        
        a_fit_pitch_R1(j,:) = a_pitch_R1;
        
        a_fit_dev_R1(j,:) = a_dev_R1;
        
        a_fit_stroke_R2(j,:) = a_stroke_R2;
        
        a_fit_pitch_R2(j,:) = a_pitch_R2;
        
        a_fit_dev_R2(j,:) = a_dev_R2;
        
        a_fit_ratio_1_2(j) = ratio_1_2;
        
        a_fit_wb1_wb2(j,1) = ceil((ratio_1_2*(n_stop(j)-n_start(j)+1))/2);
        
        a_fit_wb1_wb2(j,2) = floor(((2-ratio_1_2)*(n_stop(j)-n_start(j)+1))/2);
        
    end    
    
%         if i == 1
% 
            a_fit.stroke_L1 = a_fit_stroke_L1;
            a_fit.stroke_L2 = a_fit_stroke_L2;
            a_fit.pitch_L1 = a_fit_pitch_L1;
            a_fit.pitch_L2 = a_fit_pitch_L2;
            a_fit.dev_L1 = a_fit_dev_L1;
            a_fit.dev_L2 = a_fit_dev_L2;

            a_fit.stroke_R1 = a_fit_stroke_R1;
            a_fit.stroke_R2 = a_fit_stroke_R2;
            a_fit.pitch_R1 = a_fit_pitch_R1;
            a_fit.pitch_R2 = a_fit_pitch_R2;
            a_fit.dev_R1 = a_fit_dev_R1;
            a_fit.dev_R2 = a_fit_dev_R2;
            
            a_fit.ratio_1_2 = a_fit_ratio_1_2;
            
            a_fit.wb1_wb2 = a_fit_wb1_wb2;
            
%         else
% 
%             a_fit.stroke_L1 = [a_fit.stroke_L1; a_fit_stroke_L1];
%             a_fit.stroke_L2 = [a_fit.stroke_L2; a_fit_stroke_L2];
%             a_fit.pitch_L1 = [a_fit.pitch_L1; a_fit_pitch_L1];
%             a_fit.pitch_L2 = [a_fit.pitch_L2; a_fit_pitch_L2];
%             a_fit.dev_L1 = [a_fit.dev_L1; a_fit_dev_L1];
%             a_fit.dev_L2 = [a_fit.dev_L2; a_fit_dev_L2];
% 
%             a_fit.stroke_R1 = [a_fit.stroke_R1; a_fit_stroke_R1];
%             a_fit.stroke_R2 = [a_fit.stroke_R2; a_fit_stroke_R2];
%             a_fit.pitch_R1 = [a_fit.pitch_R1; a_fit_pitch_R1];
%             a_fit.pitch_R2 = [a_fit.pitch_R2; a_fit_pitch_R2];
%             a_fit.dev_R1 = [a_fit.dev_R1; a_fit_dev_R1];
%             a_fit.dev_R2 = [a_fit.dev_R2; a_fit_dev_R2];
%             
%             a_fit.ratio_1_2 = [a_fit.ratio_1_2; a_fit_ratio_1_2];
%             
%             a_fit.wb1_wb2 = [a_fit.wb1_wb2; a_fit_wb1_wb2];
%         end
    
    
    % Adjust wb_loc_mean and stroke_L etc. then it should work.
        
    wb_loc_mean12 = zeros(2*length(a_fit.wb1_wb2(:,1)),2);
    
    for k = 1:(length(a_fit.wb1_wb2(:,1)))
               
        wb_loc_mean12(2*k-1,:) = [1 a_fit.wb1_wb2(k,1)];
        
        wb_loc_mean12(2*k,:) = [a_fit.wb1_wb2(k,1)+1 (a_fit.wb1_wb2(k,1)+a_fit.wb1_wb2(k,2))];

    end
        
        [a_avg_stroke_L1, a_avg_stroke_L2] = average_fit(a_fit.stroke_L1,a_fit.stroke_L2,n_pol_stroke);
        
        [a_avg_stroke_R1, a_avg_stroke_R2] = average_fit(a_fit.stroke_R1,a_fit.stroke_R2,n_pol_stroke);
        
        [a_avg_pitch_L1, a_avg_pitch_L2] = average_fit(a_fit.pitch_L1,a_fit.pitch_L2,n_pol_pitch);
        
        [a_avg_pitch_R1, a_avg_pitch_R2] = average_fit(a_fit.pitch_R1,a_fit.pitch_R2,n_pol_pitch);
        
        [a_avg_dev_L1, a_avg_dev_L2] = average_fit(a_fit.dev_L1,a_fit.dev_L2,n_pol_dev);
        
        [a_avg_dev_R1, a_avg_dev_R2] = average_fit(a_fit.dev_R1,a_fit.dev_R2,n_pol_dev);        
        

        wb_loc_LR = [wb_loc_mean12; wb_loc_mean12];


        [a_avg_stroke_LR1, a_avg_stroke_LR2] = average_fit([a_fit.stroke_L1; a_fit.stroke_R1],[a_fit.stroke_L2; a_fit.stroke_R2],n_pol_stroke);

        [a_avg_pitch_LR1, a_avg_pitch_LR2] = average_fit([a_fit.pitch_L1; a_fit.pitch_R1],[a_fit.pitch_L2; a_fit.pitch_R2],n_pol_pitch);

        [a_avg_dev_LR1, a_avg_dev_LR2] = average_fit([a_fit.dev_L1; a_fit.dev_R1],[a_fit.dev_L2; a_fit.dev_R2],n_pol_dev);   
        
        
        
        % Store average coefficients in a_avg
        
        a_avg = {};

        a_avg.stroke_L1 = a_avg_stroke_L1;
        a_avg.stroke_L2 = a_avg_stroke_L2;
        a_avg.pitch_L1 = a_avg_pitch_L1;
        a_avg.pitch_L2 = a_avg_pitch_L2;
        a_avg.dev_L1 = a_avg_dev_L1;
        a_avg.dev_L2 = a_avg_dev_L2;
        

        a_avg.stroke_R1 = a_avg_stroke_R1;
        a_avg.stroke_R2 = a_avg_stroke_R2;
        a_avg.pitch_R1 = a_avg_pitch_R1;
        a_avg.pitch_R2 = a_avg_pitch_R2;
        a_avg.dev_R1 = a_avg_dev_R1;
        a_avg.dev_R2 = a_avg_dev_R2;        
        
        
        a_avg.stroke_LR1 = a_avg_stroke_LR1;
        a_avg.stroke_LR2 = a_avg_stroke_LR2;
        a_avg.pitch_LR1 = a_avg_pitch_LR1;
        a_avg.pitch_LR2 = a_avg_pitch_LR2;
        a_avg.dev_LR1 = a_avg_dev_LR1;
        a_avg.dev_LR2 = a_avg_dev_LR2;  
        
        a_avg.wb1_wb2 = [mean(a_fit.wb1_wb2(:,1)) mean(a_fit.wb1_wb2(:,2))];

    
    
    
    % Plot results:
    
    PN_stroke = Legendre_polynomial(n_pol_stroke,2,-1:0.01:1);
    PN_pitch = Legendre_polynomial(n_pol_pitch,2,-1:0.01:1);
    PN_dev = Legendre_polynomial(n_pol_dev,2,-1:0.01:1);

    nr_wb = length(a_fit.wb1_wb2(:,1))
    
    figure()
    hold on
    for k = 1:nr_wb
        m1 =a_fit.wb1_wb2(k,1);
        m2 = a_fit.wb1_wb2(k,2);
        t1 = -1:((m1/m2)/200):(-1+(m1/m2));
        t2 = (-1+(m1/m2)):((2-(m1/m2))/200):1;
        plot(t1,(PN_stroke(:,:,1)'*a_fit.stroke_L1(k,:)'),'Color',[0.5 0.5 0.5])
        plot(t2,(PN_stroke(:,:,1)'*a_fit.stroke_L2(k,:)'),'Color',[0.5 0.5 0.5])
    end
    m1 =a_avg.wb1_wb2(1);
    m2 = a_avg.wb1_wb2(2);
    t1 = -1:((m1/m2)/200):(-1+(m1/m2));
    t2 = (-1+(m1/m2)):((2-(m1/m2))/200):1;
    plot(t1,(PN_stroke(:,:,1)'*a_avg.stroke_L1),'r')
    plot(t2,(PN_stroke(:,:,1)'*a_avg.stroke_L2),'r')    
    hold off
    
    figure()
    hold on
    for k = 1:nr_wb
        m1 =a_fit.wb1_wb2(k,1);
        m2 = a_fit.wb1_wb2(k,2);
        t1 = -1:((m1/m2)/200):(-1+(m1/m2));
        t2 = (-1+(m1/m2)):((2-(m1/m2))/200):1;
        plot(t1,(PN_pitch(:,:,1)'*a_fit.pitch_L1(k,:)'),'Color',[0.5 0.5 0.5])
        plot(t2,(PN_pitch(:,:,1)'*a_fit.pitch_L2(k,:)'),'Color',[0.5 0.5 0.5])
    end
    m1 =a_avg.wb1_wb2(1);
    m2 = a_avg.wb1_wb2(2);
    t1 = -1:((m1/m2)/200):(-1+(m1/m2));
    t2 = (-1+(m1/m2)):((2-(m1/m2))/200):1;
    plot(t1,(PN_pitch(:,:,1)'*a_avg.pitch_L1),'r')
    plot(t2,(PN_pitch(:,:,1)'*a_avg.pitch_L2),'r')   
    hold off
    
    figure()
    hold on
    for k = 1:nr_wb
        m1 =a_fit.wb1_wb2(k,1);
        m2 = a_fit.wb1_wb2(k,2);
        t1 = -1:((m1/m2)/200):(-1+(m1/m2));
        t2 = (-1+(m1/m2)):((2-(m1/m2))/200):1;
        plot(-1:0.005:0,(PN_dev(:,:,1)'*a_fit.dev_L1(k,:)'),'Color',[0.5 0.5 0.5])
        plot(0:0.005:1,(PN_dev(:,:,1)'*a_fit.dev_L2(k,:)'),'Color',[0.5 0.5 0.5])
    end
    m1 =a_avg.wb1_wb2(1);
    m2 = a_avg.wb1_wb2(2);
    t1 = -1:((m1/m2)/200):(-1+(m1/m2));
    t2 = (-1+(m1/m2)):((2-(m1/m2))/200):1;
    plot(t1,(PN_dev(:,:,1)'*a_avg.dev_L1),'r')
    plot(t2,(PN_dev(:,:,1)'*a_avg.dev_L2),'r')   
    hold off
    

