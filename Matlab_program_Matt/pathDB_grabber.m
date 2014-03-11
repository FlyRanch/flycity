function output = pathDB_grabber(input)

global seqs

output = nan(1000,seqs);

for i = 1:length(input(1,:))
  k = 0;
    for j = 1:(length(input(:,1))-1)
            A = input(j,i);
            B = input(j+1,i);  
                if isfinite(A) == 1 && isfinite(B) == 1
                    k = k+1;
                    output(k,i) = A;
                end
    end
end