function [ stroke_plane ] = calc_stroke_plane_angle( tracks )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    q_left = [tracks.wingQL1,tracks.wingQL2,tracks.wingQL3,tracks.wingQL4];
    q_right = [tracks.wingQR1,tracks.wingQR2,tracks.wingQR3,tracks.wingQR4];
    %q_left = [tracks.wingQL4,tracks.wingQL1,tracks.wingQL2,tracks.wingQL3];
    %q_right = [tracks.wingQR4,tracks.wingQR1,tracks.wingQR2,tracks.wingQR3];
    track_len = size(q_left);
    track_len = track_len(1);
    left_pos(1:track_len,1:3) = nan;
    right_pos(1:track_len,1:3) = nan;
    ltip = [0,1,0];
    rtip = [0,-1,0];
    for i = 1:track_len
        left_pos(i,:) = quat2matNEW(q_left(i,:))*ltip';
        right_pos(i,:) = quat2matNEW(q_right(i,:))*rtip';
    end
    x = [right_pos(1:800,1) left_pos(1:800,1)];
    y = [right_pos(1:800,3) left_pos(1:800,3)];
    %figure();
    %plot(x,y);
    slope = Theil_Sen_Regress(x(~isnan(x)),y(~isnan(y)));
    %hold on;
    %plot(x,x*slope);
    stroke_plane = rad2deg(atan(slope))
end

