%% Turning

stroke_wb_L_MATT_bins_turning = nan(size(stroke_wb_L_MATT_bins));
stroke_wb_R_MATT_bins_turning = nan(size(stroke_wb_L_MATT_bins));
dev_wb_R_MATT_bins_turning = nan(size(stroke_wb_L_MATT_bins));
dev_wb_L_MATT_bins_turning = nan(size(stroke_wb_L_MATT_bins));
pitch_wb_R_MATT_bins_turning = nan(size(stroke_wb_L_MATT_bins));
pitch_wb_L_MATT_bins_turning = nan(size(stroke_wb_L_MATT_bins));

% stroke_wb_L_MATT_bins_turning_right = nan(size(stroke_wb_L_MATT_bins));
% stroke_wb_R_MATT_bins_turning_right = nan(size(stroke_wb_L_MATT_bins));
% dev_wb_R_MATT_bins_turning_right = nan(size(stroke_wb_L_MATT_bins));
% dev_wb_L_MATT_bins_turning_right = nan(size(stroke_wb_L_MATT_bins));
% pitch_wb_R_MATT_bins_turning_right = nan(size(stroke_wb_L_MATT_bins));
% pitch_wb_L_MATT_bins_turning_right = nan(size(stroke_wb_L_MATT_bins));
% 
% stroke_wb_L_MATT_bins_turning_left = nan(size(stroke_wb_L_MATT_bins));
% stroke_wb_R_MATT_bins_turning_left = nan(size(stroke_wb_L_MATT_bins));
% dev_wb_R_MATT_bins_turning_left = nan(size(stroke_wb_L_MATT_bins));
% dev_wb_L_MATT_bins_turning_left = nan(size(stroke_wb_L_MATT_bins));
% pitch_wb_R_MATT_bins_turning_left = nan(size(stroke_wb_L_MATT_bins));
% pitch_wb_L_MATT_bins_turning_left = nan(size(stroke_wb_L_MATT_bins));

stroke_wb_L_MATT_bins_steady = nan(size(stroke_wb_L_MATT_bins));
stroke_wb_R_MATT_bins_steady = nan(size(stroke_wb_L_MATT_bins));
dev_wb_R_MATT_bins_steady = nan(size(stroke_wb_L_MATT_bins));
dev_wb_L_MATT_bins_steady = nan(size(stroke_wb_L_MATT_bins));
pitch_wb_R_MATT_bins_steady = nan(size(stroke_wb_L_MATT_bins));
pitch_wb_L_MATT_bins_steady = nan(size(stroke_wb_L_MATT_bins));

for i = 1:length(stroke_wb_L_MATT_bins)
    
    turn = 0;
%     left_turn = 0;
%     right_turn = 0;
    steady = 0;
  
    wing = abs(stroke_wb_L_MATT_bins(:,i)-stroke_wb_R_MATT_bins(:,i));
    wing_max = max(wing);
    
     
            if wing_max > turning_tolerance
%                 left_turn = left_turn+1;
                turn = turn+1;
            else
                steady = steady+1;
            end
%             elseif wing_max < (-1*turning_tolerance)
%                 right_turn = right_turn+1;
%             elseif wing_max < turning_tolerance && wing_max > (-1*turning_tolerance)
% %                 steady = steady+1;
%             end

        if turn > steady
            
            stroke_wb_L_MATT_bins_turning(:,i) = stroke_wb_L_MATT_bins(:,i);
            stroke_wb_R_MATT_bins_turning(:,i) = stroke_wb_R_MATT_bins(:,i);
            dev_wb_R_MATT_bins_turning(:,i) = dev_wb_R_MATT_bins(:,i);
            dev_wb_L_MATT_bins_turning(:,i) = dev_wb_L_MATT_bins(:,i);
            pitch_wb_R_MATT_bins_turning(:,i) = pitch_wb_R_MATT_bins(:,i);
            pitch_wb_L_MATT_bins_turning(:,i) = pitch_wb_L_MATT_bins(:,i);
            
%     if left_turn > right_turn && left_turn > steady
%         
%         stroke_wb_L_MATT_bins_turning_left(:,i) = stroke_wb_L_MATT_bins(:,i);
%         stroke_wb_R_MATT_bins_turning_left(:,i) = stroke_wb_R_MATT_bins(:,i);
%         dev_wb_R_MATT_bins_turning_left(:,i) = dev_wb_R_MATT_bins(:,i);
%         dev_wb_L_MATT_bins_turning_left(:,i) = dev_wb_L_MATT_bins(:,i);
%         pitch_wb_R_MATT_bins_turning_left(:,i) = pitch_wb_R_MATT_bins(:,i);
%         pitch_wb_L_MATT_bins_turning_left(:,i) = pitch_wb_L_MATT_bins(:,i);
%     
%     elseif right_turn > left_turn && right_turn > steady
%         
%         stroke_wb_L_MATT_bins_turning_right(:,i) = stroke_wb_L_MATT_bins(:,i);
%         stroke_wb_R_MATT_bins_turning_right(:,i) = stroke_wb_R_MATT_bins(:,i);
%         dev_wb_R_MATT_bins_turning_right(:,i) = dev_wb_R_MATT_bins(:,i);
%         dev_wb_L_MATT_bins_turning_right(:,i) = dev_wb_L_MATT_bins(:,i);
%         pitch_wb_R_MATT_bins_turning_right(:,i) = pitch_wb_R_MATT_bins(:,i);
%         pitch_wb_L_MATT_bins_turning_right(:,i) = pitch_wb_L_MATT_bins(:,i);
        
        else
            stroke_wb_L_MATT_bins_steady(:,i) = stroke_wb_L_MATT_bins(:,i);
            stroke_wb_R_MATT_bins_steady(:,i) = stroke_wb_R_MATT_bins(:,i);
            dev_wb_R_MATT_bins_steady(:,i) = dev_wb_R_MATT_bins(:,i);
            dev_wb_L_MATT_bins_steady(:,i) = dev_wb_L_MATT_bins(:,i);
            pitch_wb_R_MATT_bins_steady(:,i) = pitch_wb_R_MATT_bins(:,i);
            pitch_wb_L_MATT_bins_steady(:,i) = pitch_wb_L_MATT_bins(:,i);
    end
end