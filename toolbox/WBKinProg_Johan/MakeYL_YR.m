function [ YL, YR ] = MakeYL_YR( qL1, qL2, qL3, qL4, qR1, qR2, qR3, qR4, start, stop, omega_on_off,dt)

    N = stop-start+1;
    
    
   if omega_on_off == 0 

       YL = zeros(4,N);

           for k = 1:N-1

               conj = [qL1(k) qL2(k) qL3(k) qL4(k)]*[qL1(k+1); qL2(k+1); qL3(k+1); qL4(k+1)];

               if conj < -0.5 && k > 1
                   qL1(k+1) = -qL1(k+1);
                   qL2(k+1) = -qL2(k+1);
                   qL3(k+1) = -qL3(k+1);
                   qL4(k+1) = -qL4(k+1);
               elseif conj < -0.5 && k == 1
                   qL1(k+1) = -qL1(k+1);
                   qL2(k+1) = -qL2(k+1);
                   qL3(k+1) = -qL3(k+1);
                   qL4(k+1) = -qL4(k+1);

                   conj2 = [qL1(k) qL2(k) qL3(k) qL4(k)]*[qL1(k+1); qL2(k+1); qL3(k+1); qL4(k+1)];

                   if conj2 < -0.5
                    qL1(k) = -qL1(k);
                    qL2(k) = -qL2(k);
                    qL3(k) = -qL3(k);
                    qL4(k) = -qL4(k);
                    
                    qL1(k+1) = -qL1(k+1);
                    qL2(k+1) = -qL2(k+1);
                    qL3(k+1) = -qL3(k+1);
                    qL4(k+1) = -qL4(k+1);
                    
                    qL1(k) = qL1(k)./norm([qL1(k); qL2(k); qL3(k); qL4(k)]);
                    qL2(k) = qL2(k)./norm([qL1(k); qL2(k); qL3(k); qL4(k)]);
                    qL3(k) = qL3(k)./norm([qL1(k); qL2(k); qL3(k); qL4(k)]);
                    qL4(k) = qL4(k)./norm([qL1(k); qL2(k); qL3(k); qL4(k)]);
                   end
               end
               
                qL1(k+1) = qL1(k+1)./norm([qL1(k+1); qL2(k+1); qL3(k+1); qL4(k+1)]);
                qL2(k+1) = qL2(k+1)./norm([qL1(k+1); qL2(k+1); qL3(k+1); qL4(k+1)]);
                qL3(k+1) = qL3(k+1)./norm([qL1(k+1); qL2(k+1); qL3(k+1); qL4(k+1)]);
                qL4(k+1) = qL4(k+1)./norm([qL1(k+1); qL2(k+1); qL3(k+1); qL4(k+1)]);
               
           end 
           


       YL(1,:) = qL1;
       YL(2,:) = qL2;
       YL(3,:) = qL3;
       YL(4,:) = qL4;
       
       
        
   
   elseif omega_on_off == 1
       
       YL = zeros(7,N);
   
       omegaL = zeros(3,N);
       
       
           for k = 2:N-1
    
                q_1 = [qL1(k-1); qL2(k-1); qL3(k-1); qL4(k-1)];

                q_1_inv = [-q_1(1); -q_1(2); -q_1(3); q_1(4)]./norm(q_1);

                q_1_inv = q_1_inv./norm(q_1_inv);

                q_2 = [qL1(k+1); qL2(k+1); qL3(k+1); qL4(k+1)];

                Q_2 = [q_2(4) q_2(3) -q_2(2) q_2(1);
                       -q_2(3) q_2(4) q_2(1) q_2(2);
                       q_2(2) -q_2(1) q_2(4) q_2(3);
                       -q_2(1) -q_2(2) -q_2(3) q_2(4)];


                q_t = Q_2*q_1_inv;

                q_t = q_t./norm(q_t);
                
                omega_norm = (2/dt)*acos(q_t(4));

                omegaL(:,k) = 0.5.*[q_t(1)*omega_norm/sin(dt*omega_norm/2); q_t(2)*omega_norm/sin(dt*omega_norm/2); q_t(3)*omega_norm/sin(dt*omega_norm/2)];
                
                clear q_1 q_1_inv q_2 Q_2 q_t R omega_norm

           end
           
          
          omegaL(:,1) = omegaL(:,2); 
          omegaL(:,N) = omegaL(:,N-1);   
       
       YL(1,:) = omegaL(1,:);
       YL(2,:) = omegaL(2,:);
       YL(3,:) = omegaL(3,:);
       YL(4,:) = qL1(:);
       YL(5,:) = qL2(:);
       YL(6,:) = qL3(:);
       YL(7,:) = qL4(:);

       
   end
   
   if omega_on_off == 0
   
       YR = zeros(4,N);

       for k = 1:N-1

           conj = [qR1(k) qR2(k) qR3(k) qR4(k)]*[qR1(k+1); qR2(k+1); qR3(k+1); qR4(k+1)];

           if conj < -0.5 && k > 1
               
               qR1(k+1) = -qR1(k+1);
               qR2(k+1) = -qR2(k+1);
               qR3(k+1) = -qR3(k+1);
               qR4(k+1) = -qR4(k+1);
               
           elseif conj < -0.5 && k == 1
               
               qR1(k+1) = -qR1(k+1);
               qR2(k+1) = -qR2(k+1);
               qR3(k+1) = -qR3(k+1);
               qR4(k+1) = -qR4(k+1);

               conj2 = [qR1(k) qR2(k) qR3(k) qR4(k)]*[qR1(k+1); qR2(k+1); qR3(k+1); qR4(k+1)];

               if conj2 < -0.5
                   
                qR1(k) = -qR1(k);
                qR2(k) = -qR2(k);
                qR3(k) = -qR3(k);
                qR4(k) = -qR4(k);
                
                qR1(k+1) = -qR1(k+1);
                qR2(k+1) = -qR2(k+1);
                qR3(k+1) = -qR3(k+1);
                qR4(k+1) = -qR4(k+1);
                
                qR1(k) = qR1(k)./norm([qR1(k); qR2(k); qR3(k); qR4(k)]);
                qR2(k) = qR2(k)./norm([qR1(k); qR2(k); qR3(k); qR4(k)]);
                qR3(k) = qR3(k)./norm([qR1(k); qR2(k); qR3(k); qR4(k)]);
                qR4(k) = qR4(k)./norm([qR1(k); qR2(k); qR3(k); qR4(k)]);
                
               end
           end    
           
        qR1(k+1) = qR1(k+1)./norm([qR1(k+1); qR2(k+1); qR3(k+1); qR4(k+1)]);
        qR2(k+1) = qR2(k+1)./norm([qR1(k+1); qR2(k+1); qR3(k+1); qR4(k+1)]);
        qR3(k+1) = qR3(k+1)./norm([qR1(k+1); qR2(k+1); qR3(k+1); qR4(k+1)]);
        qR4(k+1) = qR4(k+1)./norm([qR1(k+1); qR2(k+1); qR3(k+1); qR4(k+1)]);
       end 


       YR(1,:) = qR1;
       YR(2,:) = qR2;
       YR(3,:) = qR3;
       YR(4,:) = qR4;
   
   elseif omega_on_off == 1
       
           YR = zeros(7,N);
   
           omegaR = zeros(3,N);
           
           for k = 2:N-1
    
                q_1 = [qR1(k-1); qR2(k-1); qR3(k-1); qR4(k-1)];

                q_1_inv = [-q_1(1); -q_1(2); -q_1(3); q_1(4)]./norm(q_1);

                q_1_inv = q_1_inv./norm(q_1_inv);

                q_2 = [qR1(k+1); qR2(k+1); qR3(k+1); qR4(k+1)];

                Q_2 = [q_2(4) q_2(3) -q_2(2) q_2(1);
                       -q_2(3) q_2(4) q_2(1) q_2(2);
                       q_2(2) -q_2(1) q_2(4) q_2(3);
                       -q_2(1) -q_2(2) -q_2(3) q_2(4)];
                   
                q_t = Q_2*q_1_inv;

                q_t = q_t./norm(q_t);

                omega_norm = (2/dt)*acos(q_t(4));

                omegaR(:,k) = 0.5.*[q_t(1)*omega_norm/sin(dt*omega_norm/2); q_t(2)*omega_norm/sin(dt*omega_norm/2); q_t(3)*omega_norm/sin(dt*omega_norm/2)];
                
                clear q_1 q_1_inv q_2 Q_2 q_t R omega_norm

           end           

          omegaR(:,1) = omegaR(:,2);  
          omegaR(:,N) = omegaR(:,N-1);
       
          YR(1,:) = omegaR(1,:);
          YR(2,:) = omegaR(2,:);
          YR(3,:) = omegaR(3,:);
          YR(4,:) = qR1(:);
          YR(5,:) = qR2(:);
          YR(6,:) = qR3(:);
          YR(7,:) = qR4(:);

       
   end


end

