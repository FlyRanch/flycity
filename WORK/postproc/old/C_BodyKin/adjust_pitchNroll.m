% adjust pitch & roll
pathDB.pitch_obs = - pathDB.pitch_obs;
pathDB.pitch = - pathDB.pitch;

roll_temp = pathDB.roll_obs;
pathDB.roll_obs(roll_temp<0) = roll_temp(roll_temp<0) +180;
pathDB.roll_obs(roll_temp>0) = roll_temp(roll_temp>0) -180;
clear roll_temp

roll = pathDB.roll_obs;
t = pathDB.t;
for i=1:size(roll,2)
    roll_now = roll(:,i);

    pp = csaps(t,roll_now,.99999);
    roll_now(isnan(roll_now)==0) = fnval(pp,t(isnan(roll_now)==0));
    roll(:,i) = roll_now;    
end

pathDB.roll = roll;
clear pp i roll roll_now t

savefile = 'flightpathDB.mat';
save(savefile,'pathDB','settings')        