function output = circ_mean_deg_nonan(alpha)

alpha = alpha(isnan(alpha)==0);
alpha = deg2rad(alpha);

[mu ul ll] = circ_mean(alpha);

mu = rad2deg(mu);
ul = rad2deg(ul);
ll = rad2deg(ll);

output(1,1) = mu;
output(1,2) = ul;
output(1,3) = ll;
