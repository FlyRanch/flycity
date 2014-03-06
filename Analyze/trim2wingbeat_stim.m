function output = trim2wingbeat_stim(input)

global wingbeats_def_stim_start seqs

% output = NaN(50,50,seqs);

for i = 1:seqs
    
    count = 0;
    
    for j = 1:length(wingbeats_def_stim_start(1,:,i))
        
        if wingbeats_def_stim_start(1,j,i) > 0 && wingbeats_def_stim_start(3,j,i) < 5588
            count = count+1;
            start = wingbeats_def_stim_start(1,j,i);
            stop = wingbeats_def_stim_start(3,j,i);
            range = abs(stop-start);
            
            output(1:(range+1),count,i) = input(start:stop,i);
        end
    end
end

output(output == 0) = nan;

clear count start stop range