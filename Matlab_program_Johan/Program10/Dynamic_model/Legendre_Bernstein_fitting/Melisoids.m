function [ Mn ] = Melisoids( loc, x)

    % Create local shape functions
    
    m = length(x);
    
    n = length(loc);
    
    Mn = zeros(m, n);
    
    
    loc_max_min = zeros(n,1);
    
    for j = 1:n
        
        loc_max_min(j) = round(loc(j)*(m-1)+1);
        
    end

            
    for i = 1:n
    
        if i == 1
            
            a = loc_max_min(i):loc_max_min(i+1);
            
            delta_x = x(loc_max_min(i+1))-x(loc_max_min(i));
            
            Mn(a,i) = cos((pi/2)*((x(a)-x(loc_max_min(i)))/delta_x)).^2;
            
        elseif i == n
            
            a = loc_max_min(i-1):loc_max_min(i);
            
            delta_x = x(loc_max_min(i))-x(loc_max_min(i-1));
            
            Mn(a,i) = sin((pi/2)*((x(a)-x(loc_max_min(i-1)))/delta_x)).^2;
            
        else
            
            a = loc_max_min(i-1):loc_max_min(i);
            
            delta_x1 = x(loc_max_min(i))-x(loc_max_min(i-1));
            
            Mn(a,i) = sin((pi/2)*((x(a)-x(loc_max_min(i-1)))/delta_x1)).^2;
            
            b = (loc_max_min(i)+1):loc_max_min(i+1);
            
            delta_x2 = x(loc_max_min(i+1))-x(loc_max_min(i));
            
            Mn(b,i) = cos((pi/2)*((x(b)-x(loc_max_min(i)))/delta_x2)).^2;
            
        end
        
    end
    
end


