

rot_axis_mirror = atan2d(pitch_dot,roll_dot_mirror);

% no neg angles
rot_axis_plot = rot_axis_mirror;
rot_axis_plot(rot_axis_plot<-90) = rot_axis_plot(rot_axis_plot<-90)+360;

% remove jumps from plots (+/-180deg)
rot_axis_temp = rot_axis_plot;
for i=1:size(rot_axis_plot,2)
    for j=2:size(rot_axis_plot,1)
        if abs(rot_axis_temp(j,i) - rot_axis_temp(j-1,i)) > 80
            rot_axis_plot(j-1,i) = nan;
        end
    end
end


