function [ phi_mean, theta_mean, xsi_mean ] = Euler_mean(q1, q2, q3, q4)

    % Find the average of a quaternion by using qeneralized quaternion
    % interpolation:
    
    
    N = length(q1);
    
    % Initial estimate:
    
    m = [q1(1); q2(1); q3(1); q4(1)];
    
    err = 1;
    
    while abs(err)<0.001
        
        err_t = ones(N,1);
        
        for i = 1:N
            
            m_inv = [-m(1); -m(2); -m(3); m(4)]./(m(1)^2+m(2)^2+m(3)^2+m(4)^2);
            
            m_inv = m_inv./norm(m_inv);
            
            M_inv = [m_inv(4) m_inv(3) -m_inv(2) m_inv(1); ...
                     -m_inv(3) m_inv(4) m_inv(1) m_inv(2); ...
                     m_inv(2) -m_inv(1) m_inv(4) m_inv(3); ...
                     -m_inv(1) -m_inv(2) -m_inv(3) m_inv(4)];

            M_q = M_inv*[q1(i); q2(i); q3(i); q4(i)];
            
            M_q = M_q./norm(M_q);
                 
            err_t(i) = log(M_q);
            
            clear m_inv M_inv M_q
            
        end
        
        err = sum(err_t);
        
        m = m*exp(err);
        
        m = m./norm(m);
        
        clear err_t
        
    end
    
    phi_mean = atan2(2*(m(4)*m(1)+m(2)*m(3)),(1-2*(m(1)^2+m(2)^2)));
    theta_mean = asin(2*(m(4)*m(2)-m(3)*m(1)));
    xsi_mean = atan2(2*(m(4)*m(3)+m(1)*m(2)),1-2*(m(2)^2+m(3)^2));
    
    
    


end

