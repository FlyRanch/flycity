% Kalmar filter 3D flight track
% State = (x y z xdot ydot zdot xdotdot ydotdot zdotdot). We only observe (x y z).

% based on code that was used to generate Figure 17.9 of "Artificial Intelligence: a Modern Approach",
% Russell and Norvig, 2nd edition, Prentice Hall, in preparation.

% X(t+1) = F X(t) + noise(Q)
% Y(t) = H X(t) + noise(R)
% 

function [q,omega,settings] = batch_kalman_qEKF(q_obs,t,settings,plotting,saving)

%% batch loop
dt = t(2)-t(1);

q1_obs = q_obs(:,:,1);
q2_obs = q_obs(:,:,2);
q3_obs = q_obs(:,:,3);
q4_obs = q_obs(:,:,4);

q1 = nan(size(q1_obs));
q2 = nan(size(q1_obs));
q3 = nan(size(q1_obs));
q4 = nan(size(q1_obs));

omega_1 = nan(size(q1_obs));
omega_2 = nan(size(q1_obs));
omega_3 = nan(size(q1_obs));

%seq loop
for i=1:size(q1_obs,2)

    % measured values
    start = find(isnan(q1_obs(:,i))==0, 1 );
    stop = find(isnan(q1_obs(:,i))==0, 1, 'last' );

    counter = size(q1_obs,2) - i

    q1_raw = q1_obs(start:stop,i);
    q2_raw = q2_obs(start:stop,i);
    q3_raw = q3_obs(start:stop,i);
    q4_raw = q4_obs(start:stop,i);

    [qsmooth1,qsmooth2,qsmooth3,qsmooth4,omega_now_1,omega_now_2,omega_now_3] = filter_body(q1_raw,q2_raw,q3_raw,q4_raw,dt);
    
    q1(start:stop,i) = qsmooth1(:);
    q2(start:stop,i) = qsmooth2(:);
    q3(start:stop,i) = qsmooth3(:);
    q4(start:stop,i) = qsmooth4(:);

    omega_1(start:stop,i) = omega_now_1(:);
    omega_2(start:stop,i) = omega_now_2(:);
    omega_3(start:stop,i) = omega_now_3(:);
    
    if plotting == 1
        figure(1)

        subplot(2,1,1)
        hold off
        plot(t_now,qsmooth1,'b')
        hold on
        plot(t_now,qsmooth2,'r')
        plot(t_now,qsmooth3,'g')
        plot(t_now,qsmooth4,'c')
        grid on
        legend('q1', 'q2', 'q3', 'q4','Location','NorthWest')
        xlabel('t')
        ylabel('q')

        subplot(2,1,2)
        hold off
        plot(t_now,omega_now_1./dt,'b')
        hold on
        plot(t_now,omega_now_2./dt,'r')
        plot(t_now,omega_now_3./dt,'g')
        grid on
        legend('omega_now_1','omega_now_2','omega_now_3','Location','NorthWest')
        xlabel('t')
        ylabel('omega')

if saving == 1
        title(['kalman results ' num2str(settings.seq(i,1)) 'S' num2str(settings.seq(i,2))])
    %     saveas(gca,['xua ' num2str(settings.seq(i,1)) 'S' num2str(settings.seq(i,2)) '.fig'])
        saveas(gca,['q_' num2str(settings.seq(i,1)) 'S' num2str(settings.seq(i,2)) '.jpg'])
end
    %     pause
    end
end

q = nan(size(q_obs));
omega = nan(size(q_obs,1),size(q_obs,2),3);

q(:,:,1) = q1;
q(:,:,2) = q2;
q(:,:,3) = q3;
q(:,:,4) = q4;

omega(:,:,1) = omega_1;
omega(:,:,2) = omega_2;
omega(:,:,3) = omega_3;


