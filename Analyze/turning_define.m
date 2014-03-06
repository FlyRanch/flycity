%%% turning_def:
%
% seqs - wingbeat location - theta_R(min) - phi_R(min) - theta_L(min) -
% phi_L(min) - distance between point
%
%%%

if exist('theta_R','var') == 0
    theta_R = pathDB_grabber(pathDB.theta_R);
end

if exist('theta_L','var') == 0
    theta_L = pathDB_grabber(pathDB.theta_L);
end

turning_def = nan(length(wingbeats_def(2,:,1)),7,seqs);

for i = fliplr(1:seqs)
for j = 1:length(wingbeats_def(2,:,i))
    if isnan(wingbeats_def(2,j,i)) == 0 && isnan(wingbeats_def_L(2,j,i)) == 0
%             figure(i)
%             hold on
            
            % Right Wing
            x1 = theta_R(wingbeats_def(2,j,i));
            y1 = phi_R(wingbeats_def(2,j,i));
                        
            % Left Wing
            x2 = theta_L(wingbeats_def_L(2,j,i));
            y2 = phi_L(wingbeats_def_L(2,j,i));
                      
            turning_def(j,1,i) = i;
            turning_def(j,2,i) = wingbeats_loc(1,j,i);
            turning_def(j,3:4,i) = [x1,y1];
            turning_def(j,5:6,i) = [x2,y2];
            
            if y1 < y2
                wing{j,1,i} = 'left';
            elseif y1 > y2
                wing{j,1,i} = 'right';
            end
            
%             turning_def(j,7,i) = pdist([x1,y1; x2,y2],'euclidean');
            turning_def(j,7,i) = pdist([x1,y1; x2,y2],'euclidean');

%             plot(x1,y1,'*','color','blue','linewidth',2)
%             plot(x2,y2,'*','color','red','linewidth',2)
            
    else
        turning_def(j,1,i) = i;
        turning_def(j,2,i) = wingbeats_loc(1,j,i);
        turning_def(j,3:6,i) = nan;
        wing{j,1,i} = 'NaN';
    end
    
    turning_def(turning_def == 0) = nan;
end
end

clear input count x1 x2 y1 y2



