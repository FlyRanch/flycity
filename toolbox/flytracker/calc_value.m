function Xn = calc_value(X,n)

for i = 1:length(n)
    if isnan(n(i)) == 0
        Xn(i,1) = X(n(i),i);
    else
        Xn(i,1) = nan;
    end
end

