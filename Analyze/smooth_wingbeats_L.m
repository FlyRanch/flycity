%% Light wing

span = 12; % number of frams to average
width = 10; % number of frames to sample at each max

phi_L_smooth = phi_L;

for i = 1:seqs
    
    snip = phi_L_smooth(:,i);
    avg = nanmean(y_max_L(:,i));
    dev = nanstd(y_max_L(:,i));
    
    
    for j = 1:length(x_max_frame_L(:,i))
        
        if x_max_frame_L(j,i) > 0 && x_max_frame_L(j,i)+width <= 1000 
            if not(isnan(y_max_L(j,i))) == 1 
%                 if y_max_L(j,i) > (avg + dev)
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

phi_L = phi_L_smooth;

% right wing stroke
wingbeat_threshold = 25; % minimum number of frames that a wingbeat can be
x_max_frame_L = NaN(50,seqs);
x_min_frame_L = NaN(50,seqs);
y_max_L = NaN(50,seqs);
y_min_L = NaN(50,seqs);

for i = 1:seqs  
k=0;
m=0;
    % calculate min and max positions in frames (n) and magnitude (m)   
    for j=2:length(phi_L)-1

        A = phi_L(j,i);
        B = phi_L(j-1,i);
        C = phi_L(j+1,i);
        
        if A > B && A > C % then A is a maximum
            k=k+1;
                if k > 1
                    wingbeat_length_L = j - x_max_frame_L((k-1),i);
                    if wingbeat_length_L > wingbeat_threshold
                        y_max_L(k,i) = phi_L(j,i); % magnitude of max
                        x_max_frame_L(k,i) = (j); % frame location of min
%                         x_max_time_L(k,i) = (j)*fr; % time point of max
                    elseif wingbeat_length_L <= wingbeat_threshold
                        k=k-1; % reset k so that a zero is not left in that spot      
                    end
                    
                elseif k == 1
                        y_max_L(k,i) = phi_L(j,i); % magnitude of max
                        x_max_frame_L(k,i) = (j); % frame location of min
%                         x_max_time_L(k,i) = (j)*fr; % time point of max
                end
                
        elseif A < B && A < C % then A is a minimum
            m=m+1;
                if m > 1
                    wingbeat_length_L = j - x_min_frame_L((m-1),i);
                    if wingbeat_length_L > wingbeat_threshold        
                        y_min_L(m,i) = phi_L(j,i); % magnitude of min
                        x_min_frame_L(m,i) = (j); % frame location of min
%                         x_min_time_L(m,i) = (j)*fr; % time point of min
                    elseif wingbeat_length_L <= wingbeat_threshold
                        m = m-1; % reset m so that a zero is not left in that spot
                    end
                    
                elseif m == 1
                        y_min_L(m,i) = phi_L(j,i); % magnitude of max
                        x_min_frame_L(m,i) = (j); % frame location of min
%                         x_min_time_L(m,i) = (j)*fr; % time point of max
                end
        end
    end
end

clear A B C k m i j


%% Left wing

% Grab stroke angle data

phi_L = nan(1000,seqs);

for i = 1:length(pathDB.phi_L(1,:))
  k = 0;
    for j = 1:(length(pathDB.phi_L(:,1))-1)
            A = pathDB.phi_L(j,i);
            B = pathDB.phi_L(j+1,i);  
                if isfinite(A) == 1 && isfinite(B) == 1
                    k = k+1;
                    phi_L(k,i) = A;
                end
    end 
end

% left wing stroke
wingbeat_threshold = 25; % minimum number of frames that a wingbeat can be
x_max_frame_L = NaN(50,seqs);
x_min_frame_L = NaN(50,seqs);
y_max_L = NaN(50,seqs);
y_min_L = NaN(50,seqs);

for i = 1:seqs  
k=0;
m=0;
    % calculate min and max positions in frames (n) and magnitude (m)   
    for j=2:length(phi_L)-1

        A = phi_L(j,i);
        B = phi_L(j-1,i);
        C = phi_L(j+1,i);
        
        if A > B && A > C % then A is a maximum
            k=k+1;
                if k > 1
                    wingbeat_length_L = j - x_max_frame_L((k-1),i);
                    if wingbeat_length_L > wingbeat_threshold
                        y_max_L(k,i) = phi_L(j,i); % magnitude of max
                        x_max_frame_L(k,i) = (j); % frame location of min
%                         x_max_time_L(k,i) = (j)*fr; % time point of max
                    elseif wingbeat_length_L <= wingbeat_threshold
                        k=k-1; % reset k so that a zero is not left in that spot      
                    end
                    
                elseif k == 1
                        y_max_L(k,i) = phi_L(j,i); % magnitude of max
                        x_max_frame_L(k,i) = (j); % frame location of min
%                         x_max_time_L(k,i) = (j)*fr; % time point of max
                end
                
        elseif A < B && A < C % then A is a minimum
            m=m+1;
                if m > 1
                    wingbeat_length_L = j - x_min_frame_L((m-1),i);
                    if wingbeat_length_L > wingbeat_threshold        
                        y_min_L(m,i) = phi_L(j,i); % magnitude of min
                        x_min_frame_L(m,i) = (j); % frame location of min
%                         x_min_time_L(m,i) = (j)*fr; % time point of min
                    elseif wingbeat_length_L <= wingbeat_threshold
                        m = m-1; % reset m so that a zero is not left in that spot
                    end
                    
                elseif m == 1
                        y_min_L(m,i) = phi_L(j,i); % magnitude of max
                        x_min_frame_L(m,i) = (j); % frame location of min
%                         x_min_time_L(m,i) = (j)*fr; % time point of max
                end
        end
    end
end

clear A B C k m i j

span = 12; % number of frams to average
width = 10; % number of frames to sample at each max

phi_L_smooth = phi_L;

for i = 1:seqs
    
    snip = phi_L_smooth(:,i);
    avg = nanmean(y_max_L(:,i));
    dev = nanstd(y_max_L(:,i));
    
    
    for j = 1:length(x_max_frame_L(:,i))
        
        if x_max_frame_L(j,i) > 0 && x_max_frame_L(j,i)+width <= 1000 
            if not(isnan(y_max_L(j,i))) == 1 
%                 if y_max_L(j,i) > (avg + dev)
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

phi_L = phi_L_smooth;


%% Compare smoothed Data to original

% for i = fliplr(smooth_peaks)
% 
% figure(i)
% hold on
% title(settings.sequence_names(i))
% plot(phi_L(:,i),'-','color','blue');
% plot(phi_L_smooth(:,i),'-r');
% % plot(x_min_frame_L(:,i),y_min_L(:,i),'*','color','g')
% % plot(x_max_frame_L(:,i),y_max_L(:,i),'.','color','r')
% 
% hold off
% end

clear phi_L_smooth phi_L_smooth dev snip avg i j

