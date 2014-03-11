qbody_obs(:,:) = pathDB.qbody_obs;
qbody(:,:) = pathDB.qbody;

t = pathDB.t;


for i=1:size(qbody_obs,1)
    
    [yaw_now,pitch_now,roll_now] = quat2angle([qbody_obs(i,4) qbody_obs(i,1:3)]);
    yaw_obs(i,1) = rad2deg(yaw_now);
    pitch_obs(i,1) = rad2deg(pitch_now);
    roll_obs(i,1) = rad2deg(roll_now);

    [yaw_now,pitch_now,roll_now] = quat2angle([qbody(i,4) qbody(i,1:3)]);
    yaw(i,1) = rad2deg(yaw_now);
    pitch(i,1) = rad2deg(pitch_now);
    roll(i,1) = rad2deg(roll_now);
end

% adjust pitch & roll
pitch_obs = - pitch_obs;
roll_temp = roll_obs;
roll_obs(roll_temp<0) = roll_temp(roll_temp<0) +180;
roll_obs(roll_temp>0) = roll_temp(roll_temp>0) -180;

pitch = - pitch;
roll_temp = roll;
roll(roll_temp<0) = roll_temp(roll_temp<0) +180;
roll(roll_temp>0) = roll_temp(roll_temp>0) -180;


t = t(isnan(roll_temp)==0);

roll_obs = roll_obs(isnan(roll_temp)==0);
pitch_obs = pitch_obs(isnan(roll_temp)==0);
yaw_obs = yaw_obs(isnan(roll_temp)==0);

roll = roll(isnan(roll_temp)==0);
pitch = pitch(isnan(roll_temp)==0);
yaw = yaw(isnan(roll_temp)==0);


plot(t,roll,'r')
hold on
plot(t,pitch,'g')
plot(t,yaw,'y')
legend('roll','pitch','yaw','Location','SouthEast')

plot(t,roll_obs,'k')
plot(t,pitch_obs,'k')
plot(t,yaw_obs,'k')

plot(t,roll,'r')
plot(t,pitch,'g')
plot(t,yaw,'y')

axis tight
xlabel('t [s]')
ylabel('angle [deg]')


