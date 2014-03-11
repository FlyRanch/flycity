function [a_fit, a_avg] = WBfitting_updown(n_pol,t_loc_down,var_loc_down,t_loc_up,var_loc_up,Rds)


a_fit = {};

        [ a1, a2 ] = local_fit_updown(var_loc_down,var_loc_up,Rds,n_pol);
        
        a_fit1(j,:) = a1;
        a_fit2(j,:) = a2;
        a_fit_ratio_1_2(j) = ratio_1_2;
        
        
        a_fit_wb1_wb2(j,1) = ceil((ratio_1_2*(n_stop(j)-n_start(j)+1))/2);
        a_fit_wb1_wb2(j,2) = floor(((2-ratio_1_2)*(n_stop(j)-n_start(j)+1))/2);
    end    
    
            a_fit.var1 = a_fit1;
            a_fit.var2 = a_fit2;
            a_fit.ratio_1_2 = a_fit_ratio_1_2;
            a_fit.wb1_wb2 = a_fit_wb1_wb2;
            
    
    
    % Adjust wb_loc_mean and var then it should work.
        
    wb_loc_mean12 = zeros(2*length(a_fit.wb1_wb2(:,1)),2);
    
    for k = 1:(length(a_fit.wb1_wb2(:,1)))
               
        wb_loc_mean12(2*k-1,:) = [1 a_fit.wb1_wb2(k,1)];
        
        wb_loc_mean12(2*k,:) = [a_fit.wb1_wb2(k,1)+1 (a_fit.wb1_wb2(k,1)+a_fit.wb1_wb2(k,2))];

    end
        
        [a_avg1, a_avg2] = average_fit_periodic(a_fit.var1,a_fit.var2,n_pol);
        
        % Store average coefficients in a_avg
        
        a_avg = {};
        a_avg.var1 = a_avg1';
        a_avg.var2 = a_avg2';
        a_avg.wb1_wb2 = [mean(a_fit.wb1_wb2(:,1)) mean(a_fit.wb1_wb2(:,2))];

    
    
    
%     % Plot results:
%     
%     PN = Legendre_polynomial(n_pol,2,-1:0.01:1);
% 
%     nr_wb = length(a_fit.wb1_wb2(:,1))
%     
% %     figure()
% %     hold on
% %     for k = 1:nr_wb
% %         m1 =a_fit.wb1_wb2(k,1);
% %         m2 = a_fit.wb1_wb2(k,2);
% %         t1 = -1:((m1/m2)/200):(-1+(m1/m2));
% %         t2 = (-1+(m1/m2)):((2-(m1/m2))/200):1;
% %         plot(t1,(PN(:,:,1)'*a_fit.var1(k,:)'),'Color',[0.5 0.5 0.5])
% %         plot(t2,(PN(:,:,1)'*a_fit.var2(k,:)'),'Color',[0.5 0.5 0.5])
% %     end
%     m1 =a_avg.wb1_wb2(1);
%     m2 = a_avg.wb1_wb2(2);
%     t1 = -1:((m1/m2)/200):(-1+(m1/m2));
%     t2 = (-1+(m1/m2)):((2-(m1/m2))/200):1;
% %     plot(t1,(PN(:,:,1)'*a_avg.var1),'r')
% %     plot(t2,(PN(:,:,1)'*a_avg.var2),'r')    
%     t1_mod=(t1+1)/2;
%     t2_mod=(t2+1)/2;
%     plot(t1_mod,(PN(:,:,1)'*a_avg.var1),'r','linewidth',2)
%     plot(t2_mod,(PN(:,:,1)'*a_avg.var2),'r','linewidth',2)    
% %     hold off

