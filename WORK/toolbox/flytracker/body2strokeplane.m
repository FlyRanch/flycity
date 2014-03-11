function [xsp ysp zsp] = body2strokeplane(xb,yb,zb,strokeplane_angle)

angy = deg2rad(strokeplane_angle);

Ty = [cos(angy) 0 -sin(angy); 0 1 0; sin(angy) 0  cos(angy)];

X = Ty * [xb;yb;zb];
% X = inv(Ty) * [xb;yb;zb];

xsp = X(1);
ysp = X(2);
zsp = X(3);


