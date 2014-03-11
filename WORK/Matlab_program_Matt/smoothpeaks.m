function output = smoothpeaks(input,smooth_peaks,wingbeats_def)


output = nan(size(input));

span = 12; % number of frams to average
width = 10; % number of frames to sample at each max

for i = smooth_peaks
    snip = input(:,i);
    
    max_frame = [wingbeats_def(1,:,i), wingbeats_def(3,:,i)];
    max_frame = max_frame(isfinite(max_frame(:)));
        
    y_max = input(max_frame,i);
    
    avg = nanmean(y_max);
    dev = nanstd(y_max);
    
    for j = 1:length(max_frame)
        
        if max_frame(j) > 0 && max_frame(j)+width <= 1000
            if not(isnan(y_max(j))) == 1 
%                 if y_max_R(j,i) > (avg + dev)
                    minus = max_frame(j)-width;
                    plus = max_frame(j)+width;
                        if minus > 0
                            smoothed = smooth(snip(minus:plus),span);
                            output(minus:plus,i) = smoothed;
                        end
%                 end
            end
        end
    end
end
end