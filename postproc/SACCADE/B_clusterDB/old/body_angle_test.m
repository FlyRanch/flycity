x1 = x{1}(:);
y1 = y{1}(:);
z1 = z{1}(:);

[x2,y2,z2] = xformq_surf(x1,y1,z1,[0 0 0 0],qbody,1);
X = quatrotate(qbody,[x1,y1,z1]);

x3 = X(:,1);
y3 = X(:,2);
z3 = X(:,3);

plot3(x1,y1,z1,'.b')
hold on
plot3(x2,y2,z2,'.r')
plot3(x3,y3,z3,'.g')

axis equal
grid on
xlabel('x')
ylabel('y')
zlabel('z')

pitch = atand( (head.Position(3) -tail.Position(3)) / sqrt( (head.Position(1) -tail.Position(1))^2 + (head.Position(2) -tail.Position(2))^2 ))

[rot(1),rot(2),rot(3)] = quat2angle(qbody);
rot = rad2deg(rot)

yaw_vec = [cosd(rot(3)) sind(rot(3)) 0];
plot3([0 yaw_vec(1)],[0 yaw_vec(2)],[0 yaw_vec(3)],'k-*')

