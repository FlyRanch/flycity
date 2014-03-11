count = 0;

for i = 1:length(trims)
    A = trims(i,1);
    
    if A == 1;
        count = count+1;
        smooth_peaks(count) = i;
    end
end

smooth_peaks = fliplr(smooth_peaks);

clear count A i