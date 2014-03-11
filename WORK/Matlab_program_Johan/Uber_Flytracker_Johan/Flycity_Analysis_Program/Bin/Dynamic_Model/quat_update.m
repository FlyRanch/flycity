function [ q_new, R_mat ] = quat_update( q, w, dt )

    v_t = dt*w;
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
    q_new_t = Q_t*q;
    q_new = q_new_t/norm(q_new_t);
    
    R_mat = quat2mat(q_new);
%     R_mat = quat2mat(q_new)';

end

