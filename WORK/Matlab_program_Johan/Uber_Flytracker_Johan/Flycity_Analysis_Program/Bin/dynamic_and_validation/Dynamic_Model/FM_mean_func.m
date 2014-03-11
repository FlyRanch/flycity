function FM_mean = FM_mean_func(T,t_func,FM)

%     N = length(t_func);
% 
%     FM_func = zeros(3,N);
%     
%     for i = 1:N
% 
%         FM_func(:,i) = [interp1(T,FM(1,:),t_func(i),'spline'); ...
%                         interp1(T,FM(2,:),t_func(i),'spline'); ...
%                         interp1(T,FM(3,:),t_func(i),'spline')];
%                     
%     end
    
%     FM_func = [interp1(T,FM(1,:),t_func,'spline'); ...
%                interp1(T,FM(2,:),t_func,'spline'); ...
%                interp1(T,FM(3,:),t_func,'spline')];

    FM_func = [interp1(T,FM(1,:),t_func); ...
               interp1(T,FM(2,:),t_func); ...
               interp1(T,FM(3,:),t_func)];
    
    FM_mean = mean(FM_func,2);

end