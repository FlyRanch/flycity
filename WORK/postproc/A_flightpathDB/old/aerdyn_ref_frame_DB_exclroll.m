function [alfa, beta] = aerdyn_ref_frame_DB_exclroll(Q_body, U)

for i = 1:size(U,1)
    for j = 1:size(U,2)
        
        U_now(:,:) = U(i,j,:);
        Qbody_now(:,:) = Q_body(i,j,:);

        % Determine Euler angles of body orientation
        [yaw pitch roll] = quat2angle([Qbody_now(4);Qbody_now(1:3)]');

        % Roll transformation matrix
        T_roll = [1 0 0; 0 cos(roll) sin(roll); 0 -sin(roll) cos(roll)];

        % Remove roll transformation matrix from DCM of q_body
%         T_matrix = transpose(quat2matNEW(Qbody_now));
        T_matrix = transpose(T_roll)*transpose(quat2matNEW(Qbody_now));

        % Transform U into body reference frame
        U_body = T_matrix*U_now;

        % Determine angles alfa and beta between body reference frame and
        % aerodynamic reference frame
        alfa(i,j) = radtodeg(atan2(U_body(3),sqrt(U_body(1)^2+U_body(2)^2)));
        beta(i,j) = radtodeg(atan2(U_body(2),U_body(1)));
    end
end