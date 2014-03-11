function [spnx spny spnz] = strokeplane_normal2world(strokeplane_angle,qb)

% strokeplane2body
spn = [0;0;-1];
angy = deg2rad(strokeplane_angle);
Ty = [cos(angy) 0 -sin(angy); 0 1 0; sin(angy) 0  cos(angy)];
spnb = Ty' * spn;

% body2world
R = quat2matNEW(qb);
spnw = R * spnb;

spnx = spnw(1);
spny = spnw(2);
spnz = spnw(3);


