% [yaw, pitch, roll] = quat2angle([xh(7);xh(4:6)]');
% Y = quat2matNEW([xh(4:7)]');

[yaw, pitch, roll] = quat2angle([SOLN(1,7) SOLN(1,4:6)]);
Y = quat2matNEW([SOLN(1,4:7)]);

yaw = rad2deg(yaw)
pitch = rad2deg(pitch)
roll= 180+rad2deg(roll) 

transl = SOLN(1,1:3)'

v = [1 0 0]'
v2 = [0 -1 0]'

vz = [0 0 1]'
vy = [0 1 0]'
vx = [1 0 0]'

vy = vy*cosd(yaw) + cross(vz,vy)*sind(yaw) + (vz.'*vy)*(1-cosd(yaw))*vz
vx = vx*cosd(yaw) + cross(vz,vx)*sind(yaw) + (vz.'*vx)*(1-cosd(yaw))*vz
vx = vx*cosd(pitch) + cross(vy,vx)*sind(pitch) + (vy.'*vx)*(1-cosd(pitch))*vy

v_yaw = v*cosd(yaw) + cross(vz,v)*sind(yaw) + (vz.'*v)*(1-cosd(yaw))*vz
v2_yaw = v2*cosd(yaw) + cross(vz,v2)*sind(yaw) + (vz.'*v2)*(1-cosd(yaw))*vz

v_pitch = v_yaw*cosd(pitch) + cross(vy,v_yaw)*sind(pitch) + (vy.'*v_yaw)*(1-cosd(pitch))*vy
v2_pitch = v2_yaw*cosd(pitch) + cross(vy,v2_yaw)*sind(pitch) + (vy.'*v2_yaw)*(1-cosd(pitch))*vy

v2_roll = v2_pitch*cosd(roll) + cross(vx,v2_pitch)*sind(roll) + (vx.'*v2_pitch)*(1-cosd(roll))*vx



quiver3(transl(1),transl(2),transl(3),v_yaw(1),v_yaw(2),v_yaw(3),5,'g')
quiver3(transl(1),transl(2),transl(3),v2_yaw(1),v2_yaw(2),v2_yaw(3),5,'g')

quiver3(transl(1),transl(2),transl(3),v_pitch(1),v_pitch(2),v_pitch(3),5,'c')
quiver3(transl(1),transl(2),transl(3),v2_pitch(1),v2_pitch(2),v2_pitch(3),5,'c')

quiver3(transl(1),transl(2),transl(3),v2_roll(1),v2_roll(2),v2_roll(3),5,'k')

