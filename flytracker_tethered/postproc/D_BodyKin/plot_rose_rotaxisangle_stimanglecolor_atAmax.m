% left turn
color_var = round(abs(stim_angle_yaw_mirror_pre));

pitch_post = calc_circ_mean_value(dpitch,n_pre,n_post);
roll_post = calc_circ_mean_value(droll_mirror,n_pre,n_post);

n_now = n_Amax;
pitch_post = calc_value(dpitch,n_now);
roll_post = calc_value(droll_mirror,n_now);

rot_dir = atand(pitch_post./roll_post)
angle_post = abs(rot_dir);

figure
if exist('maxHistogramValue') == 0
    maxHistogramValue = floor(max(sqrt(pitch_post.^2 + roll_post.^2)))-1;
end

polar(0, maxHistogramValue,'-k')
hold on


for i=1:length(pitch_post)
    plot([0 pitch_post(i)],[0 roll_post(i)],'color',cmap_plot(color_var(i),:),'linewidth',2)
    hold on
end
compass_zeroup(0, 0.01,'-k');

