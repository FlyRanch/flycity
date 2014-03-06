function output = num_winbeats(input)

output = 0;

for i = 1:length(input(1,:))
    if isnan(sum(input(:,i))) == 0
        output = output+1;
    end
end
