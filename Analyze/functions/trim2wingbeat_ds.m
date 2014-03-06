function output = trim2wingbeat_ds(input)

global seqs_plot wingbeats_def trims seqs

count = 0;
output = NaN(200,(seqs*50));

for i = seqs_plot
    for j = 1:length(wingbeats_def(1,:,i))
        if wingbeats_def(1,j,i) > trims(i,2) && wingbeats_def(3,j,i) < trims(i,3)
            count = count+1;
            output(1:length(input(wingbeats_def(1,j,i):wingbeats_def(2,j,i),i)),count) = input(wingbeats_def(1,j,i):wingbeats_def(2,j,i),i);
        end
    end
end