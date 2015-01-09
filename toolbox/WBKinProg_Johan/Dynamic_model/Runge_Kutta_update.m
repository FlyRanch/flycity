function [ State_update ] = Runge_Kutta_update(Pattern,Body,State,t,dt)

    % Input: State Vector at time t.
    % Output: State Vector at time t + delta t.
    
    m = length(State);
    
    State_update = zeros(length(State),1);
   
    
    % Calculate the coefficients for the RK4 scheme:
    
    k1 = zeros(m,1);
    
    k2 = zeros(m,1);
    
    k3 = zeros(m,1);
    
    k4 = zeros(m,1);
  
    State_1 = State_space(Pattern,Body,State,t,dt);
        
%     k1(7:9) = dt*[State_1(7); State_1(8); State_1(9)];
%     
%     k1(17:19) = dt*[State_1(17); State_1(18); State_1(19)];

    k1(7:9) = [State_1(7); State_1(8); State_1(9)];
    
    k1(17:19) = [State_1(17); State_1(18); State_1(19)];

    State_2 = State_space(Pattern,Body,State+0.5*dt*k1,t+0.5*dt,dt);
%         
%     k2(7:9) = dt*[State_2(7); State_2(8); State_2(9)];
%     
%     k2(17:19) = dt*[State_2(17); State_2(18); State_2(19)];
        
    k2(7:9) = [State_2(7); State_2(8); State_2(9)];
    
    k2(17:19) = [State_2(17); State_2(18); State_2(19)]; 
    
    State_3 = State_space(Pattern,Body,State+0.5*dt*k2,t+0.5*dt,dt);
%         
%     k3(7:9) = dt*[State_3(7); State_3(8); State_3(9)];
%     
%     k3(17:19) = dt*[State_3(17); State_3(18); State_3(19)];
        
    k3(7:9) = [State_3(7); State_3(8); State_3(9)];
    
    k3(17:19) = [State_3(17); State_3(18); State_3(19)];
    
    State_4 = State_space(Pattern,Body,State+dt*k3,t+dt,dt);
%         
%     k4(7:9) = dt*[State_4(7); State_4(8); State_4(9)];
%     
%     k4(17:19) = dt*[State_4(17); State_4(18); State_4(19)];
        
    k4(7:9) = [State_4(7); State_4(8); State_4(9)];
    
    k4(17:19) = [State_4(17); State_4(18); State_4(19)];    
    
    % Update the state vector:
    
    R_body = quat2matNEW([State(10); State(11); State(12); State(13)]);
    
    State_update(1:3) = State(1:3)+R_body'*(dt*State(4:6)+(dt^2/6)*(k1(7:9)+k2(7:9)+k3(7:9)));                                                          % [mm]
    
    State_update(4:6) = State(4:6)+(dt/6)*(k1(7:9)+2*k2(7:9)+2*k3(7:9)+k4(7:9));                                                                        % [mm/s]
    
    State_update(7:9) = State(7:9);                                                                                                                     % [mm/s^2]
    
    v_t = dt*State(14:16)+(dt^2/6)*(k1(17:19)+k2(17:19)+k3(17:19));
    
    q_t = [sin(norm(v_t)/2)*(v_t(1)/norm(v_t)); sin(norm(v_t)/2)*(v_t(2)/norm(v_t)); sin(norm(v_t)/2)*(v_t(3)/norm(v_t)); cos(norm(v_t)/2)];
    
    qt2 = q_t/norm(q_t);
    
    Q_t = [ qt2(4) qt2(3) -qt2(2) qt2(1); ...
            -qt2(3) qt2(4) qt2(1) qt2(2); ...
            qt2(2) -qt2(1) qt2(4) qt2(3); ...
            -qt2(1) -qt2(2) -qt2(3) qt2(4)];

%     Q_t = [ qt2(4) -qt2(3) qt2(2) qt2(1); ...
%             qt2(3) qt2(4) -qt2(1) qt2(2); ...
%             -qt2(2) qt2(1) qt2(4) qt2(3); ...
%             -qt2(1) -qt2(2) -qt2(3) qt2(4)];

    q_new = Q_t*State(10:13);
    
    State_update(10:13) = q_new/norm(q_new);
    
    State_update(14:16) = State(14:16)+(dt/6)*(k1(17:19)+2*k2(17:19)+2*k3(17:19)+k4(17:19));
    
    State_update(17:19) = State(17:19);


    end

