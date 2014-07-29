% An_hor_max = abs(responseDB.An_hor_max);
% At_hor_max = responseDB.At_hor_max;

skip = 10

for i = 1:length(t_pre)
    
    length(t_pre)-i
    if isnan(t_pre(i)) == 0 && isnan(t_post(i)) == 0
        n_pre = find(t==t_pre(i));
        n_post = find(t==t_post(i));
        
%         for j = n_pre:skip:n_post
        for j = 200:skip:n_post

            An_hor_norm = An_hor(j,i) / max(An_hor(n_pre:n_post,i));
            At_hor_norm = At_hor(j,i) / max(At_hor(n_pre:n_post,i));

%             plot(An_hor_norm,At_hor_norm,'.')
            plot(An_hor(j,i),At_hor(j,i),'.')
            hold on
        end
    end
end
