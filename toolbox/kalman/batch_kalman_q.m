% Kalmar filter 3D flight track
% State = (x y z xdot ydot zdot xdotdot ydotdot zdotdot). We only observe (x y z).

% based on code that was used to generate Figure 17.9 of "Artificial Intelligence: a Modern Approach",
% Russell and Norvig, 2nd edition, Prentice Hall, in preparation.

% X(t+1) = F X(t) + noise(Q)
% Y(t) = H X(t) + noise(R)
% 

function [q,q_dot] = batch_kalman_q(q_obs,t,settings,plotting)



%% filter settings
filter = settings.kalman_q;

Q = eye(filter.ss);
Q(1:4,:) = filter.q*Q(1:4,:);
Q(5:8,:) = filter.q_dot*Q(5:8,:);
% Q = 0.1*eye(filter.ss);
R = filter.R*eye(filter.os);

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

q1_dot = nan(size(q1_obs));
q2_dot = nan(size(q1_obs));
q3_dot = nan(size(q1_obs));
q4_dot = nan(size(q1_obs));

%seq loop
for i=1:size(q1_obs,2)
    counter = size(q1_obs,2) - i

    % measured values
    start = find(isnan(q1_obs(:,i))==0, 1 );
    stop = find(isnan(q1_obs(:,i))==0, 1, 'last' );
    
    t_now = t(start:stop);
    Y = ([q1_obs(start:stop,i) q2_obs(start:stop,i) q3_obs(start:stop,i) q4_obs(start:stop,i)])';

    % initial state
    init_q = Y(:,1);
    init_q_dot = (Y(:,11)-Y(:,1))/10;
    
    initq = [init_q;init_q_dot];
    initq_dot = eye(filter.ss); % initial state covariance...
    % initV = 10*eye(filter.ss);

%     [xfilt, Vfilt, VVfilt, loglik] = kalman_filter(Y, filter.F, filter.H, Q, R, initq, initq_dot);
    [qsmooth, q_dotsmooth] = kalman_smoother(Y, filter.F, filter.H, Q, R, initq, initq_dot);

    q1(start:stop,i) = qsmooth(1,:);
    q2(start:stop,i) = qsmooth(2,:);
    q3(start:stop,i) = qsmooth(3,:);
    q4(start:stop,i) = qsmooth(4,:);
    
    q1_dot(start:stop,i) = qsmooth(5,:)./dt;
    q2_dot(start:stop,i) = qsmooth(6,:)./dt;
    q3_dot(start:stop,i) = qsmooth(7,:)./dt;
    q4_dot(start:stop,i) = qsmooth(8,:)./dt;

    if plotting == 1
        figure(1)
        title(['kalman results ' num2str(settings.seq(i,1)) 'S' num2str(settings.seq(i,2))])
        subplot(2,1,1)
        hold off
        plot(t_now,qsmooth(1,:),'b')
        hold on
        plot(t_now,qsmooth(2,:),'r')
        plot(t_now,qsmooth(3,:),'g')
        plot(t_now,qsmooth(4,:),'c')
        grid on
        legend('q1', 'q2', 'q3', 'q4','Location','NorthWest')
        xlabel('t')
        ylabel('q')

        subplot(2,1,2)
        hold off
        plot(t_now,qsmooth(5,:)./dt,'b')
        hold on
        plot(t_now,qsmooth(6,:)./dt,'r')
        plot(t_now,qsmooth(7,:)./dt,'g')
        plot(t_now,qsmooth(8,:)./dt,'c')
        grid on
        legend('q1_dot', 'q2_dot', 'q3_dot', 'q4_dot','Location','NorthWest')
        xlabel('t')
        ylabel('q_dot')

    %     saveas(gca,['xua ' num2str(settings.seq(i,1)) 'S' num2str(settings.seq(i,2)) '.fig'])
        saveas(gca,['q_' num2str(settings.seq(i,1)) 'S' num2str(settings.seq(i,2)) '.jpg'])
    %     pause
    end
end

q = nan(size(q_obs));
q_dot = nan(size(q_obs));

q(:,:,1) = q1;
q(:,:,2) = q2;
q(:,:,3) = q3;
q(:,:,4) = q4;

q_dot(:,:,1) = q1_dot;
q_dot(:,:,2) = q2_dot;
q_dot(:,:,3) = q3_dot;
q_dot(:,:,4) = q4_dot;
