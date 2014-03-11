function output = bin_mean(input,num)

for i = 1:200
    
        output(i,1) = mean(input(i,:));
        output(i,2) = output(i,1) + (nanstd(input(i,:))/sqrt(length(input(i,1:num))))*1.96;
        output(i,3) = output(i,1) - (nanstd(input(i,:))/sqrt(length(input(i,1:num))))*1.96;
        output(i,4) = nanstd(input(i,:));
       
end