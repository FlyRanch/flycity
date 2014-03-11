%% Right wing

span = 12; % number of frams to average
width = 10; % number of frames to sample at each max

phi_R_smooth = phi_R_stim;

for i = smooth_peaks
    
    snip = phi_R_smooth(:,i);
    avg = nanmean(y_max_R(:,i));
    dev = nanstd(y_max_R(:,i));
    
    
    for j = 1:length(x_max_frame_R(:,i))
        
        if x_max_frame_R(j,i) > 0 && x_max_frame_R(j,i)+width <= 1000 
            if not(isnan(y_max_R(j,i))) == 1 
%                 if y_max_R(j,i) > (avg + dev)
                    minus = x_max_frame_R(j,i)-width;
                    plus = x_max_frame_R(j,i)+width;
                        if minus > 0
                            smoothed = smooth(snip(minus:plus),span);
                            phi_R_smooth(minus:plus,i) = smoothed;
                        end
%                 end
            end
        end
    end
end

% Determine Max and Min of stroke angle (again, with smoothed data)

phi_R_stim = phi_R_smooth;

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
    for j=2:length(phi_R_stim)-1

        A = phi_R_stim(j,i);
        B = phi_R_stim(j-1,i);
        C = phi_R_stim(j+1,i);
        
        if A > B && A > C % then A is a maximum
            k=k+1;
                if k > 1
                    wingbeat_length_R = j - x_max_frame_R((k-1),i);
                    if wingbeat_length_R > wingbeat_threshold
                        y_max_R(k,i) = phi_R_stim(j,i); % magnitude of max
                        x_max_frame_R(k,i) = (j); % frame location of min
%                         x_max_time_R(k,i) = (j)*fr; % time point of max
                    elseif wingbeat_length_R <= wingbeat_threshold
                        k=k-1; % reset k so that a zero is not left in that spot      
                    end
                    
                elseif k == 1
                        y_max_R(k,i) = phi_R_stim(j,i); % magnitude of max
                        x_max_frame_R(k,i) = (j); % frame location of min
%                         x_max_time_R(k,i) = (j)*fr; % time point of max
                end
                
        elseif A < B && A < C % then A is a minimum
            m=m+1;
                if m > 1
                    wingbeat_length_R = j - x_min_frame_R((m-1),i);
                    if wingbeat_length_R > wingbeat_threshold        
                        y_min_R(m,i) = phi_R_stim(j,i); % magnitude of min
                        x_min_frame_R(m,i) = (j); % frame location of min
%                         x_min_time_R(m,i) = (j)*fr; % time point of min
                    elseif wingbeat_length_R <= wingbeat_threshold
                        m = m-1; % reset m so that a zero is not left in that spot
                    end
                    
                elseif m == 1
                        y_min_R(m,i) = phi_R_stim(j,i); % magnitude of max
                        x_min_frame_R(m,i) = (j); % frame location of min
%                         x_min_time_R(m,i) = (j)*fr; % time point of max
                end
        end
    end
end

clear A B C k m i j


%% Left wing

% Grab stroke angle data

for i = 1:seqs
    phi_L_stim(1:length(pathDB.phi_L(trims(i,6):end,i)),i) = pathDB.phi_L(trims(i,6):end,i);
end

% for i = 1:length(pathDB.phi_L_stim(1,:))
%   k = 0;
%     for j = 1:(length(pathDB.phi_L_stim(:,1))-1)
%             A = pathDB.phi_L_stim(j,i);
%             B = pathDB.phi_L_stim(j+1,i);  
%                 if isfinite(A) == 1 && isfinite(B) == 1
%                     k = k+1;
%                     phi_L_stim(k,i) = A;
%                 end
%     end 
% end

% left wing stroke
wingbeat_threshold = 25; % minimum number of frames that a wingbeat can be
x_max_frame_L = NaN(50,seqs);
x_min_frame_L = NaN(50,seqs);
y_max_L = NaN(50,seqs);
y_min_L = NaN(50,seqs);

for i = seqs_plot  
k=0;
m=0;
    % calculate min and max positions in frames (n) and magnitude (m)   
    for j=2:length(phi_L_stim)-1

        A = phi_L_stim(j,i);
        B = phi_L_stim(j-1,i);
        C = phi_L_stim(j+1,i);
        
        if A > B && A > C % then A is a maximum
            k=k+1;
                if k > 1
                    wingbeat_length_L = j - x_max_frame_L((k-1),i);
                    if wingbeat_length_L > wingbeat_threshold
                        y_max_L(k,i) = phi_L_stim(j,i); % magnitude of max
                        x_max_frame_L(k,i) = (j); % frame location of min
%                         x_max_time_L(k,i) = (j)*fr; % time point of max
                    elseif wingbeat_length_L <= wingbeat_threshold
                        k=k-1; % reset k so that a zero is not left in that spot      
                    end
                    
                elseif k == 1
                        y_max_L(k,i) = phi_L_stim(j,i); % magnitude of max
                        x_max_frame_L(k,i) = (j); % frame location of min
%                         x_max_time_L(k,i) = (j)*fr; % time point of max
                end
                
        elseif A < B && A < C % then A is a minimum
            m=m+1;
                if m > 1
                    wingbeat_length_L = j - x_min_frame_L((m-1),i);
                    if wingbeat_length_L > wingbeat_threshold        
                        y_min_L(m,i) = phi_L_stim(j,i); % magnitude of min
                        x_min_frame_L(m,i) = (j); % frame location of min
%                         x_min_time_L(m,i) = (j)*fr; % time point of min
                    elseif wingbeat_length_L <= wingbeat_threshold
                        m = m-1; % reset m so that a zero is not left in that spot
                    end
                    
                elseif m == 1
                        y_min_L(m,i) = phi_L_stim(j,i); % magnitude of max
                        x_min_frame_L(m,i) = (j); % frame location of min
%                         x_min_time_L(m,i) = (j)*fr; % time point of max
                end
        end
    end
end

clear A B C k m i j

span = 12; % number of frams to average
width = 10; % number of frames to sample at each max

phi_L_smooth = phi_L_stim;

for i = smooth_peaks
    
    snip = phi_L_smooth(:,i);
    avg = nanmean(y_max_L(:,i));
    dev = nanstd(y_max_L(:,i));
    
    
    for j = 1:length(x_max_frame_L(:,i))
        
        if x_max_frame_L(j,i) > 0 && x_max_frame_L(j,i)+width <= 1000 
            if not(isnan(y_max_L(j,i))) == 1 
%                 if y_max_R(j,i) > (avg + dev)
                    minus = x_max_frame_L(j,i)-width;
                    plus = x_max_frame_L(j,i)+width;
                        if minus > 0
                            smoothed = smooth(snip(minus:plus),span);
                            phi_L_smooth(minus:plus,i) = smoothed;
                        end
%                 end
            end
        end
    end
end

phi_L_stim = phi_L_smooth;


%% Compare smoothed Data to original

% for i = fliplr(smooth_peaks)
% 
% figure(i)
% hold on
% title(settings.sequence_names(i))
% plot(phi_R_stim(:,i),'-','color','blue');
% plot(phi_R_smooth(:,i),'-r');
% % plot(x_min_frame_R(:,i),y_min_R(:,i),'*','color','g')
% % plot(x_max_frame_R(:,i),y_max_R(:,i),'.','color','r')
% 
% hold off
% end

clear phi_R_smooth phi_L_smooth

