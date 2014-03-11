function PHI = EKF_PHI(x,dt)


PHI_1  = [1 0 0 0 0 0 0 0 0 0];
PHI_2  = [0 1 0 0 0 0 0 0 0 0];
PHI_3  = [0 0 1 0 0 0 0 0 0 0];
PHI_4  = [dt 0 0 1 0 0 0 0 0 0];
PHI_5  = [0 dt 0 0 1 0 0 0 0 0];
PHI_6  = [0 0 dt 0 0 1 0 0 0 0];
PHI_7  = [0 0 0 x(10)*dt/2 -x(9)*dt/2 x(8)*dt/2 1 x(6)*dt/2 -x(5)*dt/2 x(4)*dt/2];
PHI_8  = [0 0 0 x(9)*dt/2 x(10)*dt/2 -x(7)*dt/2 -x(6)*dt/2 1 x(4)*dt/2 x(5)*dt/2];
PHI_9  = [0 0 0 -x(8)*dt/2 x(7)*dt/2 x(10)*dt/2 x(5)*dt/2 -x(4)*dt/2 1 x(6)*dt/2];
PHI_10 = [0 0 0 -x(7)*dt/2 -x(8)*dt/2 -x(9)*dt/2 -x(4)*dt/2 -x(5)*dt/2 -x(6)*dt/2 1];

PHI = [PHI_1; PHI_2; PHI_3; PHI_4; PHI_5; PHI_6; PHI_7; PHI_8; PHI_9; PHI_10];
