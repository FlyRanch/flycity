function Xsp_b = rollrotx(Xsp,roll_angle)

angx = deg2rad(roll_angle);

Tx = [1 0 0; 0 cos(angx) sin(angx); 0 -sin(angx) cos(angx)];

Xsp_b = Tx * Xsp;


