function output = outlier(input,tol)

%% remove outliers
output = input;

if exist('tol') == 0
    tol = 2;
end

A = reshape(input,length(input(:,1,1)),length(input(1,:,1))*length(input(1,1,:)));

for i = 1:length(A(:,1))
    B(i,1) = nanmean(A(i,:));
    B(i,2) = nanstd(A(i,:));
end

for i = 1:length(input(1,1,:))
    
    for j = 1:length((input(:,1,i)))
        
        if isnan(nanmean(input(j,:,i))) == 0      
            for k = 1:length(input(j,:,i))        
        
                if input(j,k,i) > B(j,1) + tol*B(j,2)
                    output(j,k,i) = nan;
                end
                
                if input(j,k,i) < B(j,1) - tol*B(j,2)
                    output(j,k,i) = nan;
                end
            end
        end
    end
end
%% reshape
output = reshape(output,length(output(:,1,1)),(length(output(1,:,1))*length(output(1,1,:))));

%% 

for i = 1:length(output(1,:))
    
    if sum(isnan(output(:,i))) > 5
        
        output(:,i) = nan;
    end
end
