figure

%% initial rotation
pitch_post = pitch_dot_dot_rot1_mean/freq_steady.^2;
droll_post = roll_dot_dot_rot1_mean/freq_steady.^2;

% if exist('maxHistogramValue') == 0
    maxHistogramValue = (max(sqrt(pitch_post.^2 + droll_post.^2)));
% end

polar(0, maxHistogramValue,'-k')
hold on

for i=1:length(pitch_post)
    plot([0 pitch_post(i)],[0 droll_post(i)],'color',[0 0 1],'linewidth',2)
    hold on
end
plot([0 nanmean(pitch_post)],[0 nanmean(droll_post)],'--k','linewidth',2)

%% counter rotation
pitch_post = pitch_dot_dot_rot2_mean/freq_steady.^2;
droll_post = roll_dot_dot_rot2_mean/freq_steady.^2;

for i=1:length(pitch_post)
    plot([0 pitch_post(i)],[0 droll_post(i)],'color',[1 0 0],'linewidth',2)
    hold on
end
plot([0 nanmean(pitch_post)],[0 nanmean(droll_post)],'--k','linewidth',2)


%% average rotation axis
plot([-maxHistogramValue*sind(rot_axis_angle) maxHistogramValue*sind(rot_axis_angle)],...
     [-maxHistogramValue*cosd(rot_axis_angle) maxHistogramValue*cosd(rot_axis_angle)],...
     'color','k','linewidth',2)

plot([-maxHistogramValue*sind(rot_axis_accel_rot1rot2_mean_mean) maxHistogramValue*sind(rot_axis_accel_rot1rot2_mean_mean)],...
     [-maxHistogramValue*cosd(rot_axis_accel_rot1rot2_mean_mean) maxHistogramValue*cosd(rot_axis_accel_rot1rot2_mean_mean)],...
     '--k','linewidth',2)

compass_zeroup(0, 0.01,'-k');

