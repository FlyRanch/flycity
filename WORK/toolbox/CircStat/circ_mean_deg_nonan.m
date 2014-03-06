function [mu ul ll] = circ_mean_deg_nonan(alpha)

alpha = alpha(isnan(alpha)==0);
alpha = deg2rad(alpha);

[mu ul ll] = circ_mean(alpha);

mu = rad2deg(mu);
ul = rad2deg(ul);
ll = rad2deg(ll);
