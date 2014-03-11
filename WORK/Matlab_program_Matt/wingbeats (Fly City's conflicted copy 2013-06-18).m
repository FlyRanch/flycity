%% wingbeats
% this script uses the maximums of stroke angle for each seqs to determine
% the start and stop of each wingbeat in frames and time

global seqs wingbeats_def

%% Grab stroke angle data

phi_R = nan(1000,seqs);
offset = zeros(2,seqs);

for i = 1:length(pathDB.phi_R(1,:))
  k = 0;
    for j = 1:(length(pathDB.phi_R(:,1))-1)
            A = pathDB.phi_R(j,i);
            B = pathDB.phi_R(j+1,i);  
                if isfinite(A) == 1 && isfinite(B) == 1
                    k = k+1;
                    phi_R(k,i) = A;
                end
                
                if isfinite(A) == 0 && isfinite(B) == 1
                    offset(1,i) = j+1;
                end
                
                if isfinite(A) == 1 && isfinite(B) == 0
                    if isfinite(offset(1,i)) == 1
                        offset(2,i) = j+offset(1,i);
                    end
                end
    end 
end

clear A B k

%% Determine Max and Min of stroke angle

% right wing stroke
wingbeat_threshold = 25; % minimum number of frames that a wingbeat can be
x_max_frame_R = NaN(50,seqs);
x_min_frame_R = NaN(50,seqs);
y_max_R = NaN(50,seqs);
y_min_R = NaN(50,seqs);

for i = 1:seqs  
k=0;
m=0;
    % calculate min and max positions in frames (n) and magnitude (m)   
    for j=2:length(phi_R)-1

        A = phi_R(j,i);
        B = phi_R(j-1,i);
        C = phi_R(j+1,i);
        
        if A > B && A > C % then A is a maximum
            k=k+1;
                if k > 1
                    wingbeat_length_R = j - x_max_frame_R((k-1),i);
                    if wingbeat_length_R > wingbeat_threshold
                        y_max_R(k,i) = phi_R(j,i); % magnitude of max
                        x_max_frame_R(k,i) = (j); % frame location of min
%                         x_max_time_R(k,i) = (j)*fr; % time point of max
                    elseif wingbeat_length_R <= wingbeat_threshold
                        k=k-1; % reset k so that a zero is not left in that spot      
                    end
                    
                elseif k == 1
                        y_max_R(k,i) = phi_R(j,i); % magnitude of max
                        x_max_frame_R(k,i) = (j); % frame location of min
%                         x_max_time_R(k,i) = (j)*fr; % time point of max
                end
                
        elseif A < B && A < C % then A is a minimum
            m=m+1;
                if m > 1
                    wingbeat_length_R = j - x_min_frame_R((m-1),i);
                    if wingbeat_length_R > wingbeat_threshold        
                        y_min_R(m,i) = phi_R(j,i); % magnitude of min
                        x_min_frame_R(m,i) = (j); % frame location of min
%                         x_min_time_R(m,i) = (j)*fr; % time point of min
                    elseif wingbeat_length_R <= wingbeat_threshold
                        m = m-1; % reset m so that a zero is not left in that spot
                    end
                    
                elseif m == 1
                        y_min_R(m,i) = phi_R(j,i); % magnitude of max
                        x_min_frame_R(m,i) = (j); % frame location of min
%                         x_min_time_R(m,i) = (j)*fr; % time point of max
                end
        end
    end
end

clear A B C k m i j

%% Smooth data to get rid of clapping issues at Maxima

if exist('smooth_peaks','var') == 1

smooth_wingbeats

end
    
clear span width snip minus plus

%% Remove ouliers

% for i = 1:length(y_max_R(1,:))
%     
%     average(
%     
%     low = phi_R(:,i)<(-3);
%     high = phi_R(:,i)>(2.5);
%     
%     phi_R(low,i) = NaN;
%     phi_R(high,i) = NaN;
%     
% end
% 
% clear low high

%% Get rid of minima that are not surrounded by maxima

count = 0;
for i = 1:seqs;
    for j = 1 % length(x_min_frame_R(;,i))
        
        A = x_min_frame_R(j,i);
        B = x_max_frame_R(j,i);
        C = x_max_frame_R(j+1,i);
        
        if A < B && A < C
            x_min_frame_R(j,i) = NaN;
            y_min_R(j,i) = NaN;
            count = count+1;
        end
    end
end

if count == seqs;
    x_min_frame_R(1,:) = [];
    y_min_R(1,:) = [];
    
    A = x_min_frame_R;
    B = NaN(1,seqs);
    C = y_min_R;
    
    x_min_frame_R = [A;B];
    y_min_R = [C;B];
end

for i = 1:seqs;
    for j = 1:(length(x_min_frame_R(:,i))-1)
        
        A = x_min_frame_R(j,i);
        B = x_max_frame_R(j+1,i);
        
        if isnumeric(A) == 1 && isnumeric(B) == 1 && A > B
            x_min_frame_R(j,i) = NaN;
        end
    end
end

clear A B C g h count i j

%% Make sure maximums match real values in stroke angle
% right wing
% for i = fliplr(101:112) 
% figure(i)
% hold on
% plot(phi_R(:,i))
% 
% for j = 1:length(y_max_R(:,1))
%     
%     plot(x_max_frame_R(j,i),y_max_R(j,i),'.')
%     plot(x_min_frame_R(j,i),y_min_R(j,i),'*')
% end
% 
% 
% title(settings.sequence_names(i))
% axis([0 1000 -2.5 2.5]);
% hold off
% end

%% wingbeats - Determine wingbeat start, min, and stop

% Right wing definitons
wingbeats_def = NaN(3,50,seqs);

% Check to see if the first minima is a NaN and if so delete it
for i = 1:length(x_min_frame_R(1,:))
    
    if isnan(x_min_frame_R(1,i)) == 1
        x_min_frame_R(1:end-1,i) = x_min_frame_R(2:end,i);
    end
end

% Create wingbeat start -- min -- stop definitions
for i = 1:seqs % for each individual
    
    count = 0;
    g = -1; % start time step
    h = 0; % stop time step
    f = 1;
    num_maxima = length(x_max_frame_R(:,i));    

   for j = 1:num_maxima % for each start and stop
       
       g = g+2; % beginning
       h = h+2; % end
       f = f+2; % start of next
             
       if g < length(x_max_frame_R(:,i)) && ...
          h < length(x_max_frame_R(:,i)) && ...
          f < length(x_max_frame_R(:,i)) && ...
          x_max_frame_R(g,i) > 0 && ...
          x_max_frame_R(h,i) > 0
               
                   % wingbeat 1
                   count = count+1;
                   wingbeats_def(1,count,i) = x_max_frame_R(g,i);
                   wingbeats_def(2,count,i) = x_min_frame_R(g,i);
                   wingbeats_def(3,count,i) = x_max_frame_R(h,i);
                   
                   % wingbeat 2
                   count = count+1;
                   wingbeats_def(1,count,i) = x_max_frame_R(h,i);
                   wingbeats_def(2,count,i) = x_min_frame_R(h,i);
                   wingbeats_def(3,count,i) = x_max_frame_R(f,i);
       end
       
                if g < length(x_max_frame_R(:,i)) && x_max_frame_R(g,i) == 0;
                    g=g+1;
                end
                if h < length(x_max_frame_R(:,i)) && x_max_frame_R(h,i) == 0;
                    h=h+1;
                end
                if f < length(x_max_frame_R(:,i)) && x_max_frame_R(f,i) == 0;
                    f=f+1;
                end
   end
end

clear g h f i j count num_maxima

% Get rid of partial wingbeat defintions
for i = 1:seqs
    for j = 1:length(wingbeats_def(1,:,i))
        if isnan(sum(wingbeats_def(1:3,j,i))) == 1
            wingbeats_def(1:3,j,i) = nan;
        end
    end
end
        
% %% create wingbeats def file from offsets
% 
% wingbeats_def_offset = NaN(3,50,seqs);
% 
% for i = 1:seqs
%     for j = 1:length(wingbeats_def(:,:,i))
%         
%         A = wingbeats_def(:,j,i);
%         
%         if isfinite(A(1,1)) == 1 && isfinite(offset(1,i)) == 1
%             wingbeats_def_offset(1,j,i) = wingbeats_def(1,j,i)-(offset(1,i)+trims(i,6));
%         end
%         
%         if isfinite(A(2,1)) == 1 && isfinite(offset(1,i)) == 1
%             wingbeats_def_offset(2,j,i) = wingbeats_def(2,j,i)-(offset(1,i)+trims(i,6));
%         end
%         
%         if isfinite(A(3,1)) == 1 && isfinite(offset(1,i)) == 1
%             wingbeats_def_offset(3,j,i) = wingbeats_def(3,j,i)-(offset(1,i)+trims(i,6));
%         end
%         
%     end
% end
% 
% wingbeats_def_offset - wingbeats_def;


clear smoothed wingbeat_length_L wingbeat_length_R wingbeat_min_R wingbeat_start_R wingbeat_stop_R wingbeat_threshold
clear wingbeats_def_offset x_max_frame_L x_max_frame_R x_min_frame_R x_min_frame_L y_max_R y_max_L y_min_L y_min_R



