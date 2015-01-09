
% rotation vel from roll_dot & pitch_dot
% rot_dot_L_mirror = (roll_dot_mirror - pitch_dot) * cosd(rot_axis_angle);
% rot_dot_R_mirror = (roll_dot_mirror + pitch_dot) * cosd(rot_axis_angle);
% 
rot_dot_L_mirror = roll_dot_mirror*sind(rot_axis_angle) - pitch_dot*cosd(rot_axis_angle);
rot_dot_R_mirror = roll_dot_mirror*cosd(rot_axis_angle) + pitch_dot*sind(rot_axis_angle);

rot_axis_mirror = atan2d(pitch_dot,roll_dot_mirror);

%% angular accels by num diff
rot_dot_dot_L_mirror = nan(size(rot_dot_L_mirror));
rot_dot_dot_R_mirror = nan(size(rot_dot_R_mirror));

rot_dot_dot_L_mirror(2:end-1,:) = rot_dot_L_mirror(3:end,:) - rot_dot_L_mirror(1:end-2,:);
rot_dot_dot_R_mirror(2:end-1,:) = rot_dot_R_mirror(3:end,:) - rot_dot_R_mirror(1:end-2,:);

rot_dot_dot_L_mirror = rot_dot_dot_L_mirror/2/dt;
rot_dot_dot_R_mirror = rot_dot_dot_R_mirror/2/dt;

%% angular pos by num int
drot_L_mirror = nan(size(rot_dot_L_mirror));
drot_R_mirror = nan(size(rot_dot_R_mirror));

for j = 1:size(rot_dot_L_mirror,2)
    for i = 2:size(rot_dot_L_mirror,1)
        if isnan(rot_dot_L_mirror(i,j))==0
            
                if isnan(drot_L_mirror(i-1,j))==0
                    drot_L_mirror_prev = drot_L_mirror(i-1,j);
                    drot_R_mirror_prev = drot_R_mirror(i-1,j);
                else
                    drot_L_mirror_prev = 0;
                    drot_R_mirror_prev = 0;
                end
                % integrate
                drot_L_mirror(i,j) = drot_L_mirror_prev + rot_dot_L_mirror(i,j)*dt;
                drot_R_mirror(i,j) = drot_R_mirror_prev + rot_dot_R_mirror(i,j)*dt;
        end
    end
end

drot_L_plot = drot_L_mirror;
rot_dot_L_plot = rot_dot_L_mirror;
rot_dot_dot_L_plot = rot_dot_dot_L_mirror;

drot_R_plot = drot_R_mirror;
rot_dot_R_plot = rot_dot_R_mirror;
rot_dot_dot_R_plot = rot_dot_dot_R_mirror;

% no angles < -90
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

