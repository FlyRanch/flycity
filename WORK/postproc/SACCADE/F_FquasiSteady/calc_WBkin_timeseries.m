function [stroke,stroke_dot,stroke_dot_dot,dev,dev_dot,dev_dot_dot,rot,rot_dot,rot_dot_dot,dt] = ...
    calc_WBkin_timeseries(t,stroke,dev,rot)

%  calc time vars
dt = gradient(t);

%% derivatives
stroke_dot     = gradient(stroke)./dt;
dev_dot     = gradient(dev)./dt;
rot_dot     = gradient(rot)./dt;

stroke_dot_dot     = gradient(stroke_dot)./dt;
dev_dot_dot     = gradient(dev_dot)./dt;
rot_dot_dot     = gradient(rot_dot)./dt;

%% deg2rad
stroke = deg2rad(stroke);
stroke_dot = deg2rad(stroke_dot);
stroke_dot_dot = deg2rad(stroke_dot_dot);

dev = deg2rad(dev);
dev_dot = deg2rad(dev_dot);
dev_dot_dot = deg2rad(dev_dot_dot);

rot = deg2rad(rot);
rot_dot = deg2rad(rot_dot);
rot_dot_dot = deg2rad(rot_dot_dot);
