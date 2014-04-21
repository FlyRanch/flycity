An_hor_max = abs(responseDB.An_hor_max);
At_hor_max = responseDB.At_hor_max;

An_hor_norm = nan(size(An_hor));
At_hor_norm = nan(size(At_hor));
skip = 100
for j = 1:length(An_hor_max)
    length(An_hor_max)-j
    
    An_hor_norm(:,j) = An_hor(:,j) / An_hor_max(j);
    At_hor_norm(:,j) = At_hor(:,j) / At_hor_max(j);
end

for i = 1:length(t_pre)
    if isnan(t_pre(i)) == 0 && isnan(t_post(i)) == 0
        n_pre = find(t==t_pre(i));
        n_post = find(t==t_post(i));
        
        for k = n_pre:skip:n_post
            plot(An_hor_norm(k,i),At_hor_norm(k,i),'.')
            hold on
        end
    end
end
