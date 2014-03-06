%% numbers wingbeats with respect to stimulus

wingbeats_loc = wingbeats_def(1,:,:);

for i = 1:length(wingbeats_def(1,:,:))
    wingbeats_loc(1,:,i) = (wingbeats_def(1,:,i)+offset(1,i));
end   

stimulus = trims(:,6); % location of the stimulus frame

% average_wingbeat_length = nanmean(nanmean(reshape((wingbeats_def(3,:,:)-wingbeats_def(1,:,:)),length(wingbeats_def(1,:,1)),length(wingbeats_def(1,1,:)))));

for i = 1:seqs
    
    count = 0;
    
    stim_location = trims(i,6); % location of stimulus
    values = (wingbeats_loc(1,:,i)-stim_location); % wingbeats locations in reference to stim
    values(values < 0) = nan; % remove wingbeats before stimulus
    
    [I,J] = nanmin(values);
    
    for j = J:50
        count = count+1;
        wingbeats_loc(1,j,i) = count;
    end
    
    I = 1:(J-1);
    
    if length(I(1,:)) > 1
        
        I(2,:) = fliplr(1:max(I))*-1;
        
        for k = 1:length(I(1,:))
            wingbeats_loc(1,k,i) = I(2,k);
        end
    end   
end
        
clear values stim_location I J i j k

% reallign wingbeats to locations

minimum = min(wingbeats_loc(:,1,:));

for i = 1:length(wingbeats_loc(1,:,1))
wingbeats_locations(1,i,1:seqs) = minimum;
minimum = minimum + 1;
end

wingbeats_loc(2,:,:) = nan;
for i = 1:seqs
wingbeats_loc(2,:,i) = linspace(1,length(wingbeats_loc(2,:,i)),length(wingbeats_loc(2,:,i)));
end

wingbeats_definitions = nan(3,length(wingbeats_loc(1,:,1)),seqs);

for i = 1:seqs
    count = 0;
    
    J = min(min(wingbeats_loc(:,1,:))):max(max(wingbeats_loc(1,:,:)));
    J(J~=0);
    
    for j = J
        
        if abs(j) > 0
            count = count+1;
            if isfinite(find(wingbeats_loc(1,:,i) == j)) == 1
                column = find(wingbeats_loc(1,:,i) == j);
                            
                wingbeats_definitions(1:3,count,i) = wingbeats_def(1:3,column,i);

            end
        end
    end
end

wingbeats_definitions(:,51:end,:) = [];

wingbeats_def = wingbeats_definitions;
wingbeats_loc = wingbeats_locations;
clear wingbeats_definitions wingbeats_locations number minimum loc ans J

% wingbeats_def(:,7,:) = [];
% wingbeats_def(:,50,:) = nan;

for i = 1:seqs
wingbeats_loc(1,7:49,i) = wingbeats_loc(1,8:50,i);
wingbeats_loc(1,50,i) = wingbeats_loc(1,49,i)+1;
end