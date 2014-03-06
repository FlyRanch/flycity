U = sqrt(pathDB.vel(:,:,1).^2 + pathDB.vel(:,:,2).^2 + pathDB.vel(:,:,3).^2 );
At = pathDB.At;
An = pathDB.An;
alpha_dot = pathDB.alpha_dot;

At_hor = pathDB.At_hor;
An_hor = pathDB.An_hor;
alpha_dot_hor = pathDB.alpha_dot_hor;

pitch_incroll = pathDB.pitch_incroll;
pitch_NOroll = -pathDB.pitch_NOroll;
pitch_global = pathDB.pitch_global;
roll = pathDB.roll;
t = pathDB.t;

trigger_frame = find(t == min(abs(t)));

% steady state body angle
b=[62 53 37 28 22];
u=[0 .2 .4 .6 .9];
plot(u,b)
ub = csaps(u,b);
plot(u,b)
hold on
fnplt(ub)
b_steady = fnval(ub,U);

Dpitch_incroll = pitch_incroll - b_steady;
Dpitch_NOroll = pitch_NOroll - b_steady;
Dpitch_global = pitch_global - b_steady;

figure
hold on
plot(Dpitch_incroll(trigger_frame,:),At(trigger_frame,:),'.b')
plot(Dpitch_NOroll(trigger_frame,:),At(trigger_frame,:),'.r')
plot(Dpitch_global(trigger_frame,:),At(trigger_frame,:),'.k')

pp=polyfit(Dpitch_global(trigger_frame,:),At(trigger_frame,:),1)
x=[-20,15]
y=polyval(pp,x)
plot(x,y,'-k')

figure
hold on
plot(Dpitch_incroll(1:100:trigger_frame,:),At(1:100:trigger_frame,:),'.b')
plot(Dpitch_NOroll(1:100:trigger_frame,:),At(1:100:trigger_frame,:),'.r')
plot(Dpitch_global(1:100:trigger_frame,:),At(1:100:trigger_frame,:),'.k')

pp=polyfit(Dpitch_global(1:100:trigger_frame,:),At(1:100:trigger_frame,:),1)
x=[-20,15]
y=polyval(pp,x)
plot(x,y,'-k')

figure
hold on
plot(Dpitch_incroll(1:100:end,:),At_hor(1:100:end,:),'.b')
plot(Dpitch_NOroll(1:100:end,:),At_hor(1:100:end,:),'.r')
plot(Dpitch_global(1:100:end,:),At_hor(1:100:end,:),'.k')

pp=polyfit(Dpitch_global(1:100:trigger_frame,:),At(1:100:trigger_frame,:),1)
x=[-20,15]
y=polyval(pp,x)
plot(x,y,'-k')

figure
hold on
plot(pitch_incroll(1:100:trigger_frame,:),At_hor(1:100:trigger_frame,:),'.b')
plot(pitch_NOroll(1:100:trigger_frame,:),At_hor(1:100:trigger_frame,:),'.r')
plot(pitch_global(1:100:trigger_frame,:),At_hor(1:100:trigger_frame,:),'.k')


figure
hold on
plot(abs(roll(1:100:end,:)),An(1:100:end,:),'.r')
plot(abs(roll(1:100:end,:)),An_hor(1:100:end,:),'.k')

figure
hold on
plot(abs(roll(1:100:trigger_frame,:)),alpha_dot(1:100:trigger_frame,:),'.r')
plot(abs(roll(1:100:trigger_frame,:)),alpha_dot_hor(1:100:trigger_frame,:),'.k')

