function [ q_bar ] = quat_mean( t_func, T, q_raw )


    N = length(t_func);

    q_func = zeros(4,N);
    
    for i = 1:N

        q_t = [interp1(T,q_raw(1,:),t_func(i),'spline'); ...
               interp1(T,q_raw(2,:),t_func(i),'spline'); ...
               interp1(T,q_raw(3,:),t_func(i),'spline'); ...
               interp1(T,q_raw(4,:),t_func(i),'spline')];

        q_func(:,i) = q_t/norm(q_t);
                    
    end
    
%     [ q_bar ] = q_avg(q_func(1,:), q_func(2,:) ,q_func(3,:), q_func(4,:));

    q_bar = avg_quaternion_markley(q_func');

end

