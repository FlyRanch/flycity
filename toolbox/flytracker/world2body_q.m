function [xb yb zb] = world2body_q(xw,yw,zw,qb)

R = quat2matNEW(qb);
X = R' * [xw;yw;zw];

xb = X(1);
yb = X(2);
zb = X(3);


