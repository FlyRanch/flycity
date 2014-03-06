function output = trim(input)

output = input;


for i = 1:length(input(:,1))
    
    A = input(i,:);
    B = nanmean(A);
    C = nanstd(A);
    
    for j = 1:length(A);
        
        if input(i,j) > (B + 2*C)
            output(i,j) = nan;
        end
        
        if input(i,j) < (B - 2*C)
            output(i,j) = nan;
        end
    end
    
end
    
    
    