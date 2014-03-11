% left turn
color_var = round(abs(turn_angle_vel_mirror));

pitch_post = calc_circ_mean_value(dpitch,n_pre,n_post);
droll_post = calc_circ_mean_value(droll_mirror,n_pre,n_post);
% roll_post = calc_circ_mean_value(roll_mirror,n_pre,n_post);

% n_now = n_Amax;
% pitch_post = calc_value(dpitch,n_now);
% droll_post = calc_value(droll_mirror,n_now);

% rot_dir = atand(pitch_post./droll_post)
rot_dir = atan2d(pitch_post,droll_post)
angle_post = abs(rot_dir);

figure
if exist('maxHistogramValue') == 0
    maxHistogramValue = floor(max(sqrt(pitch_post.^2 + droll_post.^2)))-1;
end

polar(0, maxHistogramValue,'-k')
hold on


for i=1:length(pitch_post)
    plot([0 pitch_post(i)],[0 droll_post(i)],'color',cmap_plot(color_var(i),:),'linewidth',2)
    hold on
end
compass_zeroup(0, 0.01,'-k');

