function [ L_wingtip_path, R_wingtip_path, L_wingbeat_loc, R_wingbeat_loc] = Wingbeat_LocFromStroke_direct( qL1, qL2, qL3, qL4, qR1, qR2, qR3, qR4, wing_l)

    % Calculate and return: wingtip path per wingbeat and the time of
    % stroke reversal (seen sideways from the body), the first time is the start of the first measured
    % downstroke.
    
    N = length(qL1);
    
    Lwt = wing_l.*[0; -1; 0]; % location left wingtip
    Rwt = wing_l.*[0; 1; 0];  % location right wingtip
            
    % Calculate left wing_tip path
    
    for j = 1:N

        DCM_L = quat2matNEW([qL1(j) qL2(j) qL3(j) qL4(j)]);

        Lwingtip(j,:) = DCM_L*Lwt;
        
        clear DCM_L

    end
    
    % Calculate right wing_tip path
    
    for j = 1:N

        DCM_R = quat2matNEW([qR1(j) qR2(j) qR3(j) qR4(j)]);

        Rwingtip(j,:) = DCM_R*Rwt;
        
        clear DCM_R

    end
        
    % Find the locations of the start of the upstroke and the start of the
    % downstroke (sequences will start with the first downstroke)
    
    [GammaL_dwn, locL_dwn] = findpeaks(Lwingtip(:,3), 'minpeakdistance' , 20);
    [GammaL_up, locL_up] = findpeaks(-Lwingtip(:,3), 'minpeakdistance' , 20);
    
%     LocGamma = sortrows([locL_dwn GammaL_dwn; locL_up GammaL_up],1);
%     locL = LocGamma(:,1);
%     GammaL = LocGamma(:,2);
    
    [GammaR_dwn, locR_dwn] = findpeaks((Rwingtip(:,3)), 'minpeakdistance' , 20);
    [GammaR_up, locR_up] = findpeaks((-Rwingtip(:,3)), 'minpeakdistance' , 20);
    
%     LocGamma = sortrows([locR_dwn GammaR_dwn; locR_up GammaR_up],1);
%     locR = LocGamma(:,1);
%     GammaR = LocGamma(:,2);

%% old code
% 
%     if length(locL) ~= length(locR)
%         
%         'missed stroke'           
%         length(locL)
%         length(locR)
% 
%         % fix too many wb
%         wb_lengthL = locL(2:end) - locL(1:end-1);
%         i = 0;
%         while min(wb_lengthL) < 10
%             i=i+1;
%             if wb_lengthL(i) < 10
%                 locL(i) = [];
%                 GammaL(i) = [];
%                 i=0;
%                 wb_lengthL = locL(2:end) - locL(1:end-1);
%             end
%         end
% 
%         wb_lengthR = locR(2:end) - locR(1:end-1);
%         i = 0;
%         while min(wb_lengthR) < 10
%             i=i+1;
%             if wb_lengthR(i) < 10
%                 locR(i) = [];
%                 GammaR(i) = [];
%                 i=0;
%                 wb_lengthR = locR(2:end) - locR(1:end-1);
%             end
%         end
% 
%         % fix missing wb
%         LRdiff = length(locL) - length(locR);
%         
%         stop_loop = 0;
%         while LRdiff ~= 0 && stop_loop<100
%             
%             stop_loop = stop_loop+1;
% 
%             for i = 1:min([length(locL) length(locR)])
%                 loc_diff = locL(i) - locR(i);
% 
%                 if loc_diff > 5
%                     locL = [locL(1:i-1);locR(i);locL(i:end)];
%                     GammaL = [GammaL(1:i-1);GammaR(i);GammaL(i:end)];
%                 elseif loc_diff < -5
%                     locR = [locR(1:i-1);locL(i);locR(i:end)];
%                     GammaR = [GammaR(1:i-1);GammaL(i);GammaR(i:end)];
%                 end
%             end
%             
%             LRdiff = length(locL) - length(locR);
% 
%             % remove last
%             if abs(loc_diff) < 5 && LRdiff ~= 0
%                 if LRdiff > 0
%                     locL(end) = [];
%                     GammaL(end) = [];
%                 elseif LRdiff < 0
%                     locR(end) = [];
%                     GammaR(end) = [];
%                 end
%             end
%                 
%             LRdiff = length(locL) - length(locR);
%         end
%     end
%    
%     
%     
%    
%     
%     if abs(sum(GammaL(1:2:end))) > abs(sum(GammaL(2:2:end)))
%        
%         loc_upst_L = locL(1:2:end);
%         loc_dwnst_L = locL(2:2:end);
%         
%     elseif abs(sum(GammaL(1:2:end))) < abs(sum(GammaL(2:2:end)))
%        
%         loc_upst_L = locL(2:2:end);
%         loc_dwnst_L = locL(1:2:end);
%     else
%        
%         loc_upst_L = 1;
%         loc_dwnst_L = 1;
%         
%     end
%     
% %     if length(loc_upst_L) > length(loc_dwnst_L)
% %         
% %         loc_upst_L = [loc_upst_L(1:end-1)];
% %         
% %     elseif length(loc_upst_L) < length(loc_dwnst_L)
% %         
% %         loc_dwnst_L = [ loc_dwnst_L(1:end-1) ];
% %         
% %     end
% 
%         
%     if abs(sum(GammaR(1:2:end))) > abs(sum(GammaR(2:2:end)))
%        
%         loc_upst_R = locR(1:2:end);
%         loc_dwnst_R = locR(2:2:end);
%         
%     elseif abs(sum(GammaR(1:2:end))) < abs(sum(GammaR(2:2:end)))
%        
%         loc_upst_R = locR(2:2:end);
%         loc_dwnst_R = locR(1:2:end);
%         
%     else
%         loc_upst_R = 1;
%         loc_dwnst_R = 1;
%     end
%     
% %     if length(loc_upst_R) > length(loc_dwnst_R)
% %         
% %         loc_upst_R = [ loc_upst_R(1:end-1) ];
% %         
% %     elseif length(loc_upst_R) < length(loc_dwnst_R)
% %         
% %         loc_dwnst_R = [ loc_dwnst_R(1:end-1) ];
% %         
% %     end
%     
%         %Check whether all the peaks are found
%     
%     if max(abs(loc_upst_L-loc_upst_R)) > 8
%         'missed upstroke'
%         
%         
% %         id_missed = find(abs(loc_upst_L-loc_upst_R)> 8)
% %         loc_upst_L
% %         loc_upst_R
% %         abs(loc_upst_L-loc_upst_R)
%     elseif max(abs(loc_dwnst_L-loc_dwnst_R)) > 8
%         'missed downstroke'
% %         abs(loc_dwnst_L-loc_dwnst_R)
% %         loc_dwnst_L
% %         loc_dwnst_R
%     end    
% %     
% %     %Check whether all the peaks are found
% %     
% %     if max(abs(loc_upst_L-loc_upst_R)) > 8
% %         'missed upstroke'
% %         abs(loc_upst_L-loc_upst_R)
% %     elseif max(abs(loc_dwnst_L-loc_dwnst_R)) > 8
% %         'missed downstroke'
% %         abs(loc_dwnst_L-loc_dwnst_R)
% %     end
%    
%     
%     
%     
%     % Return the wingbeat location and the wingtip path:
%     
%     L_wingtip_path = Lwingtip;
%     
%     R_wingtip_path = Rwingtip;
%     
%     L_wingbeat_loc = nan(max(length(loc_upst_L),length(loc_dwnst_L)),2);
%     
%     R_wingbeat_loc = nan(max(length(loc_upst_R),length(loc_dwnst_R)),2);
%     
%     L_wingbeat_loc(1:length(loc_upst_L),1) = loc_upst_L;
%     
%     L_wingbeat_loc(1:length(loc_dwnst_L),2) = loc_dwnst_L;
%     
%     R_wingbeat_loc(1:length(loc_upst_R),1) = loc_upst_R;
%     
%     R_wingbeat_loc(1:length(loc_dwnst_R),2) = loc_dwnst_R;
    
    

    
    %% Return the wingbeat location and the wingtip path:
    
    L_wingtip_path = Lwingtip;
    
    R_wingtip_path = Rwingtip;
    
    L_wingbeat_loc = nan(max(length(locL_up),length(locL_dwn)),2);
    
    R_wingbeat_loc = nan(max(length(locR_up),length(locR_dwn)),2);
    
    L_wingbeat_loc(1:length(locL_up),1) = locL_up;
    
    L_wingbeat_loc(1:length(locL_dwn),2) = locL_dwn;
    
    R_wingbeat_loc(1:length(locR_up),1) = locR_up;
    
    R_wingbeat_loc(1:length(locR_dwn),2) = locR_dwn;
end

