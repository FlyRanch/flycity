function [xb yb zb] = world2body_euler(xw,yw,zw,angx,angy,angz)

angx = deg2rad(angx);
angy = deg2rad(angy);
angz = deg2rad(angz);

Tx = [1 0 0; 0 cos(angx) sin(angx); 0 -sin(angx) cos(angx)];
Ty = [cos(angy) 0 -sin(angy); 0 1 0; sin(angy) 0  cos(angy)];
Tz = [cos(angz) sin(angz) 0; -sin(angz) cos(angz) 0; 0 0 1];

T = Tx * Ty * Tz;

% X = T * [xw;yw;zw];
X = T' * [xw;yw;zw];

xb = X(1);
yb = X(2);
zb = X(3);


