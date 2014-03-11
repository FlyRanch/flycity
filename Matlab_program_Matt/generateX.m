function output = generateX(input)

global seqs

output = NaN(200,(seqs*50));

for i = 1:length(input)
A = sum(not(isnan(input(:,i))));
output(1:A,i) = linspace(0,1,A);
end