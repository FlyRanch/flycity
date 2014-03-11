function [Y] = MakeY(q1, q2, q3, q4, N, dt, omega_on_off)

    % First check whether the quaternion representation is smooth and
    % doesn't has jumps to the non-existing negative quaternion. The
    % following loop corrects the quaternions if necessary.
        
    if omega_on_off == 0
        
    Y = zeros(4,N);

  
   for k = 1:N-1
       
       conj = [q1(k) q2(k) q3(k) q4(k)]*[q1(k+1); q2(k+1); q3(k+1); q4(k+1)];
       
       if conj < -0.5 && k>1
           q1(k+1) = -q1(k+1);
           q2(k+1) = -q2(k+1);
           q3(k+1) = -q3(k+1);
           q4(k+1) = -q4(k+1);
       elseif conj < -0.5 && k == 1
           q1(k+1) = -q1(k+1);
           q2(k+1) = -q2(k+1);
           q3(k+1) = -q3(k+1);
           q4(k+1) = -q4(k+1);
           
           conj2 = [q1(k) q2(k) q3(k) q4(k)]*[q1(k+1); q2(k+1); q3(k+1); q4(k+1)];
           
           if conj2 < -0.5
            q1(k) = -q1(k);
            q2(k) = -q2(k);
            q3(k) = -q3(k);
            q4(k) = -q4(k);
            
            q1(k+1) = -q1(k+1);
            q2(k+1) = -q2(k+1);
            q3(k+1) = -q3(k+1);
            q4(k+1) = -q4(k+1);
            
            q1(k) = q1(k)./norm([q1(k); q2(k); q3(k); q4(k)]);
            q2(k) = q2(k)./norm([q1(k); q2(k); q3(k); q4(k)]);
            q3(k) = q3(k)./norm([q1(k); q2(k); q3(k); q4(k)]);
            q4(k) = q4(k)./norm([q1(k); q2(k); q3(k); q4(k)]);
           end
       end    
       
       q1(k+1) = q1(k+1)./norm([q1(k+1); q2(k+1); q3(k+1); q4(k+1)]);
       q2(k+1) = q2(k+1)./norm([q1(k+1); q2(k+1); q3(k+1); q4(k+1)]);
       q3(k+1) = q3(k+1)./norm([q1(k+1); q2(k+1); q3(k+1); q4(k+1)]);
       q4(k+1) = q4(k+1)./norm([q1(k+1); q2(k+1); q3(k+1); q4(k+1)]);
       
   end 
   
   Y(1,:) = q1(:);
   Y(2,:) = q2(:);
   Y(3,:) = q3(:);
   Y(4,:) = q4(:);
   
   
   
   elseif omega_on_off == 1
        
        
        
        
    Y = zeros(7,N);
   
    omega = zeros(3,N);
    
    for k = 2:N-1
    
    q_1 = [q1(k-1); q2(k-1); q3(k-1); q4(k-1)];
    
    q_1_inv = [-q_1(1); -q_1(2); -q_1(3); q_1(4)]./norm(q_1);
    
    q_1_inv = q_1_inv./norm(q_1_inv);
    
    q_2 = [q1(k+1); q2(k+1); q3(k+1); q4(k+1)];
    
    Q_2 = [q_2(4) q_2(3) -q_2(2) q_2(1);
           -q_2(3) q_2(4) q_2(1) q_2(2);
           q_2(2) -q_2(1) q_2(4) q_2(3);
           -q_2(1) -q_2(2) -q_2(3) q_2(4)];     
            
    q_t = Q_2*q_1_inv;
    
    q_t = q_t./norm(q_t);
    
    omega_norm = (2/dt)*acos(q_t(4)); 
    
    omega(:,k) = 0.5.*[q_t(1)*omega_norm/sin(dt*omega_norm/2); q_t(2)*omega_norm/sin(dt*omega_norm/2); q_t(3)*omega_norm/sin(dt*omega_norm/2)];
    
    end
    
   
    omega(:,1) = omega(:,2);
    omega(:,N) = omega(:,N-1);
   
   Y(1,:) = omega(1,:);
   Y(2,:) = omega(2,:);
   Y(3,:) = omega(3,:);
   Y(4,:) = q1(:);
   Y(5,:) = q2(:);
   Y(6,:) = q3(:);
   Y(7,:) = q4(:);
   


   
   
    end

