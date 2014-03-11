function [s s0] = circ_std_deg_nonan(alpha)

alpha = alpha(isnan(alpha)==0);
alpha = deg2rad(alpha);

[s s0] = circ_std(alpha);

s = rad2deg(s);
s0 = rad2deg(s0);

