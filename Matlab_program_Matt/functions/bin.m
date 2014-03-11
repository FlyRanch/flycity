function output = bin(input,n)

input_unwrap = unwrap(deg2rad(input)); % convert to radians and unwrap

output = nan(n,length(input(1,:)));


for i = 1:length(input_unwrap(1,:))
    
    start = 1;
    
    for j = 1:length(input_unwrap(:,i))-1 % determine start and stop of seqs
        
        A = isnan(input_unwrap(j,i));
        B = isnan(input_unwrap(j+1,i));
        
        if A == 0 && B == 1;
            stop = j;
        end
    end
    
    if exist('stop','var') == 1
        
        input_now = input_unwrap(start:stop,i); % trim to seqs length (remove nans)
        n_input_now = [1:(length(input_now)-1)/(n-1):length(input_now)]; % 
        input_now_interp = rad2deg(interp1(input_now,n_input_now));
    
        output(1:length(input_now_interp),i) = input_now_interp;
    end
end

% wraptopi(output);

end

