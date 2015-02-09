
% allign data
allign_angle = angle_vel_trig2resp;
allign_angle = calc_value(angle_vel,max(n_first,n_resp-100));
for i=1:size(pos_x,2)
    pos_x_align(:,i) =  pos_x(:,i) .* cosd(allign_angle(i)) -  pos_y(:,i) .* sind(allign_angle(i));
    pos_y_align(:,i) =  pos_x(:,i) .* sind(allign_angle(i)) +  pos_y(:,i) .* cosd(allign_angle(i));
end
pos_z_align =  pos_z;

% center data
pos_x_align_pre = calc_value(pos_x_align,n_pre);
pos_y_align_pre = calc_value(pos_y_align,n_pre);
pos_z_align_pre = calc_value(pos_z_align,n_pre);
for i=1:size(pos_x,2)
    pos_x_align(:,i) = pos_x_align(:,i) - pos_x_align_pre(i);
    pos_y_align(:,i) = pos_y_align(:,i) - pos_y_align_pre(i);
    pos_z_align(:,i) = pos_z_align(:,i) - pos_z_align_pre(i);
    
    
%     plot(pos_x_align(:,i),pos_y_align(:,i),'r')
%     hold on
%     axis equal
%     pause
    
end

% mirror data
pos_x_align_mirror = pos_x_align;
for i=1:size(pos_x,2)
    
    pos_x_align_now = pos_x_align(:,i);
    pos_x_last = pos_x_align_now(find((isnan(pos_x_align_now)==0), 1, 'last' ));
    
    if pos_x_last < 0
        pos_x_align_mirror(:,i) = -pos_x_align_mirror(:,i);
    end
end

% frame stamp
n = [1:size(V,1)]';
n_pos = [];
for i=1:size(V,2)    
    n_pos = [n_pos n-n_pre(i)];
end
n_pos = n_pos.*pos_x*0+n_pos;

cm=[0:.01:1];
cm0=zeros(size(cm));
cm1=ones(size(cm));
cmap=[cm0 cm ; cm0 cm0 ; 1-cm cm0]';
cmap=[0.5*cm 0.5*(1+cm) ; .5*cm1 .5*cm1 ; .5*(2-cm) .5*(1-cm)]';

cmap_steps=length(cm)-1;

n_start = t_start/dt;
n_stop = t_stop/dt;
n_pos_plot = n_pos;
n_pos_plot(n_pos_plot<n_start)=nan;
n_pos_plot(n_pos_plot>n_stop)=nan;

pos_on = n_pos*0+1;
pos_x_align_mirror_plot = pos_x_align_mirror.*pos_on;
pos_x_align_plot = pos_x_align.*pos_on;
pos_y_align_plot = pos_y_align.*pos_on;
pos_z_align_plot = pos_z_align.*pos_on;

%% plot paths
% plot data color:time NO mirror
if make_mov == 1

% 2D
skip_fig = 40;

close all
figure(1)
hold on
axis equal

figure(2)
hold on
axis equal

figure(3)
hold on
axis equal

n_col = n_pos_plot - min(n_pos_plot(:))+1;
col_2d = round(n_col./max(abs(n_col(:)))*2*cmap_steps);
col_2d(col_2d<1)=1;

for i=1:size(pos_x_align_plot,2)
    for j=1:skip_fig:size(pos_x_align_plot,1)-skip_fig
    
        if isnan(pos_x_align_plot(j,i)) == 0 && isnan(pos_x_align_plot(j+skip_fig,i)) == 0 && isnan(col_2d(j,i)) == 0
            
             % top view
             figure(1)
             plot([pos_x_align_plot(j,i)*1000 pos_x_align_plot(j+skip_fig,i)*1000],...
            [pos_y_align_plot(j,i)*1000 pos_y_align_plot(j+skip_fig,i)*1000],...
            '-','color',cmap(col_2d(j,i),:),'linew',2)
    
             % front view
             figure(2)
             plot([pos_x_align_plot(j,i)*1000 pos_x_align_plot(j+skip_fig,i)*1000],...
            [pos_z_align_plot(j,i)*1000 pos_z_align_plot(j+skip_fig,i)*1000],...
            '-','color',cmap(col_2d(j,i),:),'linew',2)
    
             % side view
             figure(3)
             plot([pos_y_align_plot(j,i)*1000 pos_y_align_plot(j+skip_fig,i)*1000],...
            [pos_z_align_plot(j,i)*1000 pos_z_align_plot(j+skip_fig,i)*1000],...
            '-','color',cmap(col_2d(j,i),:),'linew',2)
    
        end
    end
end

figure(1)
axis tight
colormap(cmap)
axis([-10 10 -10 10])
set(gca,'XTick',-10:5:10) 
set(gca,'YTick',-10:5:10,'fontsize',20)
xlabel('y [mm]','fontsize',20)
ylabel('x [mm]','fontsize',20)

saveas(gca,['MSfigs_flightpaths_time',num2str(t_start),'to',num2str(t_stop),'s_mirror_topview.fig'])
saveas(gca,['MSfigs_flightpaths_time',num2str(t_start),'to',num2str(t_stop),'s_mirror_topview.png'])
plot2svg(['MSfigs_flightpaths_time',num2str(t_start),'to',num2str(t_stop),'s_mirror_topview.svg'])

figure(2)
axis tight
colormap(cmap)
axis([-10 10 -5 5])
set(gca,'XTick',-10:5:10) 
set(gca,'YTick',-10:5:10,'fontsize',20) 
xlabel('x [mm]','fontsize',20)
ylabel('z [mm]','fontsize',20)

saveas(gca,['MSfigs_flightpaths_time',num2str(t_start),'to',num2str(t_stop),'s_mirror_frontview.fig'])
saveas(gca,['MSfigs_flightpaths_time',num2str(t_start),'to',num2str(t_stop),'s_mirror_frontview.png'])
plot2svg(['MSfigs_flightpaths_time',num2str(t_start),'to',num2str(t_stop),'s_mirror_frontview.svg'])

figure(3)
axis tight
colormap(cmap)
axis([-10 10 -5 5])
set(gca,'XTick',-10:5:10) 
set(gca,'YTick',-10:5:10,'fontsize',20) 
xlabel('y [mm]','fontsize',20)
ylabel('z [mm]','fontsize',20)

saveas(gca,['MSfigs_flightpaths_time',num2str(t_start),'to',num2str(t_stop),'s_mirror_sideview.fig'])
saveas(gca,['MSfigs_flightpaths_time',num2str(t_start),'to',num2str(t_stop),'s_mirror_sideview.png'])
plot2svg(['MSfigs_flightpaths_time',num2str(t_start),'to',num2str(t_stop),'s_mirror_sideview.svg'])

figure(2)
colorbar
saveas(gca,['MSfigs_flightpaths_colorbar.fig'])
saveas(gca,['MSfigs_flightpaths_colorbar.png'])
plot2svg(['MSfigs_flightpaths_colorbar.svg'])

% 3D surfaces
figure
hold on
axis equal
surface([0 0;0 0],[0 0;0 0],[0 0;0 0],[-1 -1;1 1],'facecol','no','edgecol','no')
% colorbar

col = n_pos_plot./max(abs(n_pos_plot(:)));
for i=1:size(pos_x_align_plot,2)
    
%     % 2D path
%     surface([pos_x_align_plot(:,i)*1000 pos_x_align_plot(:,i)*1000],...
%         [pos_y_align_plot(:,i)*1000 pos_y_align_plot(:,i)*1000],...
%         [zeros(size(pos_x_align_plot(:,i))) zeros(size(pos_x_align_plot(:,i)))],...
%         [col(:,i) col(:,i)],'facecol','no','edgecol','interp','linew',2)

    % 3D path
    surface([pos_x_align_plot(:,i)*1000 pos_x_align_plot(:,i)*1000],...
        [pos_y_align_plot(:,i)*1000 pos_y_align_plot(:,i)*1000],...
        [pos_z_align_plot(:,i)*1000 pos_z_align_plot(:,i)*1000],...
        [col(:,i) col(:,i)],'facecol','no','edgecol','interp','linew',2)

%     plot(pos_x_align_plot(:,i)*1000,pos_y_align_plot(:,i)*1000,'r')
%     pause
end
axis tight
colormap(cmap)
axis([-10 10 -10 10 -5 5])
set(gca,'XTick',-10:5:10) 
set(gca,'YTick',-10:5:10) 
set(gca,'ZTick',-10:5:10,'fontsize',8)

xlabel('x [mm]')
ylabel('y [mm]')
zlabel('z [mm]','fontsize',10)
grid on

% view(2)
% saveas(gca,['MSfigs_flightpaths_time',num2str(t_start),'to',num2str(t_stop),'s_topview.fig'])
% saveas(gca,['MSfigs_flightpaths_time',num2str(t_start),'to',num2str(t_stop),'s_topview.png'])
% % plot2svg(['MSfigs_flightpaths_time',num2str(t_start),'to',num2str(t_stop),'s_topview.png'])

view(3)
view(-135,30)
view(-225,30)
saveas(gca,['MSfigs_flightpaths_time',num2str(t_start),'to',num2str(t_stop),'s_persp1.fig'])
saveas(gca,['MSfigs_flightpaths_time',num2str(t_start),'to',num2str(t_stop),'s_persp1.png'])
% plot2svg(['MSfigs_flightpaths_time',num2str(t_start),'to',num2str(t_stop),'s_persp1.png'])

view(3)
view(-135,30)
% view(-225,30)
saveas(gca,['MSfigs_flightpaths_time',num2str(t_start),'to',num2str(t_stop),'s_persp2.fig'])
saveas(gca,['MSfigs_flightpaths_time',num2str(t_start),'to',num2str(t_stop),'s_persp2.png'])
% plot2svg(['MSfigs_flightpaths_time',num2str(t_start),'to',num2str(t_stop),'s_persp2.png'])

% view(90,0)
% saveas(gca,['MSfigs_flightpaths_time',num2str(t_start),'to',num2str(t_stop),'s_sideview.fig'])
% saveas(gca,['MSfigs_flightpaths_time',num2str(t_start),'to',num2str(t_stop),'s_sideview.png'])
% % plot2svg(['MSfigs_flightpaths_time',num2str(t_start),'to',num2str(t_stop),'s_sideview.png'])

% rotation movie
mkdir('saccade_flightpaths')
cd('saccade_flightpaths')
axis off
for i=1:100;
    az = 3.6*i + 135;
    view(az,30)
%     pause(.04)
    axis equal
    axis([-10 10 -10 10 -5 5])
    saveas(gca,['saccade_flightpaths_',num2str(i),'.tif'])

end
cd ..

%% plot data color:time INC mirror
% 2D
skip_fig = 40;

close all
figure(1)
hold on
axis equal

figure(2)
hold on
axis equal

figure(3)
hold on
axis equal

n_col = n_pos_plot - min(n_pos_plot(:))+1;
col_2d = round(n_col./max(abs(n_col(:)))*2*cmap_steps);
col_2d(col_2d<1)=1;

for i=1:size(pos_x_align_mirror_plot,2)
    for j=1:skip_fig:size(pos_x_align_mirror_plot,1)-skip_fig
    
        if isnan(pos_x_align_mirror_plot(j,i)) == 0 && isnan(pos_x_align_mirror_plot(j+skip_fig,i)) == 0 && isnan(col_2d(j,i)) == 0
            
             % top view
             figure(1)
%              plot([pos_x_align_mirror_plot(j,i)*1000 pos_x_align_mirror_plot(j+skip_fig,i)*1000],...
%             [pos_y_align_plot(j,i)*1000 pos_y_align_plot(j+skip_fig,i)*1000],...
%             '-','color',cmap(col_2d(j,i),:),'linew',2)
    
             plot([pos_y_align_plot(j,i)*1000 pos_y_align_plot(j+skip_fig,i)*1000],...
            [pos_x_align_mirror_plot(j,i)*1000 pos_x_align_mirror_plot(j+skip_fig,i)*1000],...
            '-','color',cmap(col_2d(j,i),:),'linew',2)
    
             % front view
             figure(2)
             plot([pos_x_align_mirror_plot(j,i)*1000 pos_x_align_mirror_plot(j+skip_fig,i)*1000],...
            [pos_z_align_plot(j,i)*1000 pos_z_align_plot(j+skip_fig,i)*1000],...
            '-','color',cmap(col_2d(j,i),:),'linew',2)
    
             % side view
             figure(3)
             plot([pos_y_align_plot(j,i)*1000 pos_y_align_plot(j+skip_fig,i)*1000],...
            [pos_z_align_plot(j,i)*1000 pos_z_align_plot(j+skip_fig,i)*1000],...
            '-','color',cmap(col_2d(j,i),:),'linew',2)
    
        end
    end
end

figure(1)
axis tight
colormap(cmap)
axis([-10 10 -1 10])
set(gca,'XTick',-10:5:10) 
set(gca,'YTick',-10:5:10,'fontsize',20)
xlabel('y [mm]','fontsize',20)
ylabel('x [mm]','fontsize',20)

saveas(gca,['MSfigs_flightpaths_time',num2str(t_start),'to',num2str(t_stop),'s_mirror_topview.fig'])
saveas(gca,['MSfigs_flightpaths_time',num2str(t_start),'to',num2str(t_stop),'s_mirror_topview.png'])
plot2svg(['MSfigs_flightpaths_time',num2str(t_start),'to',num2str(t_stop),'s_mirror_topview.svg'])

figure(2)
axis tight
colormap(cmap)
axis([-1 10 -5 5])
set(gca,'XTick',-10:5:10) 
set(gca,'YTick',-10:5:10,'fontsize',20) 
xlabel('x [mm]','fontsize',20)
ylabel('z [mm]','fontsize',20)

saveas(gca,['MSfigs_flightpaths_time',num2str(t_start),'to',num2str(t_stop),'s_mirror_frontview.fig'])
saveas(gca,['MSfigs_flightpaths_time',num2str(t_start),'to',num2str(t_stop),'s_mirror_frontview.png'])
plot2svg(['MSfigs_flightpaths_time',num2str(t_start),'to',num2str(t_stop),'s_mirror_frontview.svg'])

figure(3)
axis tight
colormap(cmap)
axis([-10 10 -5 5])
set(gca,'XTick',-10:5:10) 
set(gca,'YTick',-10:5:10,'fontsize',20) 
xlabel('y [mm]','fontsize',20)
ylabel('z [mm]','fontsize',20)

saveas(gca,['MSfigs_flightpaths_time',num2str(t_start),'to',num2str(t_stop),'s_mirror_sideview.fig'])
saveas(gca,['MSfigs_flightpaths_time',num2str(t_start),'to',num2str(t_stop),'s_mirror_sideview.png'])
plot2svg(['MSfigs_flightpaths_time',num2str(t_start),'to',num2str(t_stop),'s_mirror_sideview.svg'])

% 3D
figure
hold on
axis equal
surface([0 0;0 0],[0 0;0 0],[0 0;0 0],[-1 -1;1 1],'facecol','no','edgecol','no')
% colorbar

col = n_pos_plot./max(abs(n_pos_plot(:)));
for i=1:size(pos_x_align_mirror_plot,2)
    
    % 3D path
    surface([pos_x_align_mirror_plot(:,i)*1000 pos_x_align_mirror_plot(:,i)*1000],...
        [pos_y_align_plot(:,i)*1000 pos_y_align_plot(:,i)*1000],...
        [pos_z_align_plot(:,i)*1000 pos_z_align_plot(:,i)*1000],...
        [col(:,i) col(:,i)],'facecol','no','edgecol','interp','linew',2)

%     plot(pos_x_align_plot(:,i)*1000,pos_y_align_plot(:,i)*1000,'r')
%     pause
end
axis tight
colormap(cmap)
axis([-1 10 -10 10 -5 5])
set(gca,'XTick',-10:5:10) 
set(gca,'YTick',-10:5:10) 
set(gca,'ZTick',-10:5:10,'fontsize',8)

xlabel('x [mm]')
ylabel('y [mm]')
zlabel('z [mm]','fontsize',10)
grid on

% view(2)
% saveas(gca,['MSfigs_flightpaths_time',num2str(t_start),'to',num2str(t_stop),'s_mirror_topview.fig'])
% saveas(gca,['MSfigs_flightpaths_time',num2str(t_start),'to',num2str(t_stop),'s_mirror_topview.png'])
% % plot2svg(['MSfigs_flightpaths_time',num2str(t_start),'to',num2str(t_stop),'s_mirror_topview.png'])

view(3)
view(-135,30)
view(-225,30)
saveas(gca,['MSfigs_flightpaths_time',num2str(t_start),'to',num2str(t_stop),'s_mirror_persp1.fig'])
saveas(gca,['MSfigs_flightpaths_time',num2str(t_start),'to',num2str(t_stop),'s_mirror_persp1.png'])
% plot2svg(['MSfigs_flightpaths_time',num2str(t_start),'to',num2str(t_stop),'s_mirror_persp1.png'])

view(3)
view(-135,30)
% view(-225,30)
saveas(gca,['MSfigs_flightpaths_time',num2str(t_start),'to',num2str(t_stop),'s_mirror_persp2.fig'])
saveas(gca,['MSfigs_flightpaths_time',num2str(t_start),'to',num2str(t_stop),'s_mirror_persp2.png'])
% plot2svg(['MSfigs_flightpaths_time',num2str(t_start),'to',num2str(t_stop),'s_mirror_persp2.png'])

% view(90,0)
% saveas(gca,['MSfigs_flightpaths_time',num2str(t_start),'to',num2str(t_stop),'s_mirror_sideview.fig'])
% saveas(gca,['MSfigs_flightpaths_time',num2str(t_start),'to',num2str(t_stop),'s_mirror_sideview.png'])
% % plot2svg(['MSfigs_flightpaths_time',num2str(t_start),'to',num2str(t_stop),'s_mirror_sideview.png'])

% rotation movie
mkdir('saccade_flightpaths_mirror')
cd('saccade_flightpaths_mirror')
axis off
for i=1:100;
    az = 3.6*i + 135;
    view(az,30)
%     pause(.04)
    axis equal
    axis([-1 10 -10 10 -5 5])
    saveas(gca,['saccade_flightpaths_mirror_',num2str(i),'.tif'])

end
cd ..

%% make movies
% if make_mov == 1
    
dn_mov = 3;
mov_start = min(n_pos_plot(:));
mov_stop = max(n_pos_plot(:));

% movie no mirror accelColor
cmap_steps = 100;
Amax = 10;
% Amax = max(abs(A(:)));

cmap_mov=jet(cmap_steps);
col = round(A./Amax*cmap_steps).*pos_on;
col(col>cmap_steps)=cmap_steps;

mkdir('saccades_topview_accelColor')
cd('saccades_topview_accelColor')
% movfilename = ['saccades_topview_accelColor.avi'];
% mov = avifile(movfilename, 'compression', 'none');

figure
hold on
axis equal
axis([min(pos_x_align_plot(:)) max(pos_x_align_plot(:)) min(pos_y_align_plot(:)) max(pos_y_align_plot(:))])
axis off
for i=mov_start:dn_mov:mov_stop
    
    x_now = pos_x_align_plot(n_pos_plot==i);
    y_now = pos_y_align_plot(n_pos_plot==i);
    col_now = col(n_pos_plot==i);

    for j=1:length(x_now)    
        if isnan(col_now(j))==0
            plot(x_now(j),y_now(j),'.','color',cmap_mov(col_now(j),:))
        end
    end
    
%     if i<0
%         plot(x_now,y_now,'.b')
%     else
%         plot(x_now,y_now,'.r')
%     end
    
%     mov = addframe(mov,gcf);
%     axis equal
%     axis([-10 10 -10 10 -5 5])
    saveas(gcf,['saccades_topview_accelColor_Amax',num2str(Amax),'mps2_',num2str(i-mov_start+1),'.tif']);
%     pause(.001)
end
% mov = close(mov);
cd ..

% movie INC mirror accelColor
mkdir('saccades_topview_accelColor_mirror')
cd('saccades_topview_accelColor_mirror')
% movfilename = ['saccades_topview_accelColor_mirror.avi'];
% mov = avifile(movfilename, 'compression', 'none');

figure
hold on
axis equal
axis([min(pos_x_align_mirror_plot(:)) max(pos_x_align_mirror_plot(:)) min(pos_y_align_plot(:)) max(pos_y_align_plot(:))])
axis off

for i=mov_start:dn_mov:mov_stop
    
    x_now = pos_x_align_mirror_plot(n_pos_plot==i);
    y_now = pos_y_align_plot(n_pos_plot==i);
    col_now = col(n_pos_plot==i);

    for j=1:length(x_now)    
        if isnan(col_now(j))==0
            plot(x_now(j),y_now(j),'.','color',cmap_mov(col_now(j),:))
        end
    end

%     if i<0
%         plot(x_now,y_now,'.b')
%     else
%         plot(x_now,y_now,'.r')
%     end
    
%     mov = addframe(mov,gcf);
%     axis equal
%     axis([-1 10 -10 10 -5 5])
    saveas(gcf,['saccades_topview_accelColor_mirror_Amax',num2str(Amax),'mps2_',num2str(i-mov_start+1),'.tif']);
%     pause(.001)
end
% mov = close(mov);
cd ..

% movie no mirror VelColor
cmap_steps = 100;
Vmax = 0.3;
% Vmax = max(abs(V(:)));

cmap_mov=jet(cmap_steps);
col = round(V./Vmax*cmap_steps).*pos_on;
col(col>cmap_steps)=cmap_steps;

mkdir('saccades_topview_VelColor')
cd('saccades_topview_VelColor')
% movfilename = ['saccades_topview_VelColor.avi'];
% mov = avifile(movfilename, 'compression', 'none');

figure
hold on
axis equal
axis([min(pos_x_align_plot(:)) max(pos_x_align_plot(:)) min(pos_y_align_plot(:)) max(pos_y_align_plot(:))])
axis off
for i=mov_start:dn_mov:mov_stop
    
    x_now = pos_x_align_plot(n_pos_plot==i);
    y_now = pos_y_align_plot(n_pos_plot==i);
    col_now = col(n_pos_plot==i);

    for j=1:length(x_now)    
        if isnan(col_now(j))==0
            plot(x_now(j),y_now(j),'.','color',cmap_mov(col_now(j),:))
        end
    end
    
%     if i<0
%         plot(x_now,y_now,'.b')
%     else
%         plot(x_now,y_now,'.r')
%     end
    
%     mov = addframe(mov,gcf);
    saveas(gcf,['saccades_topview_VelColor_Vmax',num2str(Vmax),'mps_',num2str(i-mov_start+1),'.tif']);
%     pause(.001)
end
% mov = close(mov);
cd ..

% movie INC mirror VelColor
mkdir('saccades_topview_VelColor_mirror')
cd('saccades_topview_VelColor_mirror')
% movfilename = ['saccades_topview_VelColor_mirror.avi'];
% mov = avifile(movfilename, 'compression', 'none');

figure
hold on
axis equal
axis([min(pos_x_align_mirror_plot(:)) max(pos_x_align_mirror_plot(:)) min(pos_y_align_plot(:)) max(pos_y_align_plot(:))])
axis off

for i=mov_start:dn_mov:mov_stop
    
    x_now = pos_x_align_mirror_plot(n_pos_plot==i);
    y_now = pos_y_align_plot(n_pos_plot==i);
    col_now = col(n_pos_plot==i);

    for j=1:length(x_now)    
        if isnan(col_now(j))==0
            plot(x_now(j),y_now(j),'.','color',cmap_mov(col_now(j),:))
        end
    end

%     if i<0
%         plot(x_now,y_now,'.b')
%     else
%         plot(x_now,y_now,'.r')
%     end
    
%     mov = addframe(mov,gcf);
    saveas(gcf,['saccades_topview_VelColor_mirror_Vmax',num2str(Vmax),'mps_',num2str(i-mov_start+1),'.tif']);
%     pause(.001)
end
% mov = close(mov);
cd ..

end
