function output = trim2wingbeat_wb(input,seqs)

global wingbeats_def


% output = NaN(200,50,seqs));

for i = 1:seqs
    
    count = 0;
    
    for j = 1:length(wingbeats_def(1,:,i))
        count = count+1;
        if wingbeats_def(1,j,i) > 0 && wingbeats_def(3,j,i) < 5588
            
%             count = count+1;
            start = wingbeats_def(1,j,i);
            stop = wingbeats_def(3,j,i);
            range = abs(stop-start);
            
            output(1:(range+1),count,i) = input(start:stop,i);
        end
    end
end

output(output == 0) = nan;

clear count start stop range