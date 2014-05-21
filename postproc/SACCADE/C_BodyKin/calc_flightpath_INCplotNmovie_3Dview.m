
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
if plot_paths == 1
    
    mkdir('flightpaths')
    cd('flightpaths')

%% movie TimeColor

cmap_mov = cmap;
n_col = n_pos_plot - min(n_pos_plot(:))+1;
col = round(n_col./max(abs(n_col(:)))*2*cmap_steps);
cmap_mov(end/2,:)=[];
col=col+1;

figure
colormap(cmap_mov)
colorbar
saveas(gcf,['colorbar_Time.fig']);
saveas(gcf,['colorbar_Time.tif']);
plot2svg(['colorbar_Time.svg']);

close all
figure(1)
hold on
axis equal
axis([-10 10 -10 10])
axis off

figure(2)
hold on
axis equal
axis([-10 10 -5 5])
axis off

figure(3)
hold on
axis equal
axis([-10 10 -5 5])
axis off

figure(4)
hold on
axis equal
axis([-10 10 -10 10 -5 5])
axis off

figure(5)
hold on
axis equal
axis([-1 10 -10 10])
axis off

figure(6)
hold on
axis equal
axis([-1 10 -5 5])
axis off

figure(7)
hold on
axis equal
axis([-10 10 -5 5])
axis off

figure(8)
hold on
axis equal
axis([-1 10 -10 10 -5 5])
axis off

if make_mov == 1
    mkdir('saccades_TimeColor')
    cd('saccades_TimeColor')

    figure(1)
    saveas(gcf,['saccades_topview_TimeColor',num2str(0),'.tif']);

    figure(2)
    saveas(gcf,['saccades_frontview_TimeColor',num2str(0),'.tif']);

    figure(3)
    saveas(gcf,['saccades_sideview_TimeColor',num2str(0),'.tif']);

    figure(4)
    view(-225,30)
    saveas(gcf,['saccades_perspview1_TimeColor',num2str(0),'.tif']);
    view(-135,30)
    saveas(gcf,['saccades_perspview2_TimeColor',num2str(0),'.tif']);

    figure(5)
    saveas(gcf,['saccades_topview_TimeColor_mirror',num2str(0),'.tif']);

    figure(6)
    saveas(gcf,['saccades_frontview_TimeColor_mirror',num2str(0),'.tif']);

    figure(7)
    saveas(gcf,['saccades_sideview_TimeColor_mirror',num2str(0),'.tif']);

    figure(8)
    view(-225,30)
    saveas(gcf,['saccades_perspview1_TimeColor_mirror',num2str(0),'.tif']);
    view(-135,30)
    saveas(gcf,['saccades_perspview2_TimeColor_mirror',num2str(0),'.tif']);
end

if make_mov == 1
    dn_mov = 3;
else
    dn_mov = 40;
end

mov_start = min(n_pos_plot(:));
mov_stop = max(n_pos_plot(:));

for i=mov_start:dn_mov:mov_stop
    
    [n_now,m_now] = find(n_pos_plot==i);

    for j=1:length(n_now)
        if isnan(col(n_now(j),m_now(j)))==0 && isnan(col(n_now(j)+dn_mov,m_now(j)))==0
            
            figure(1)
            plot([pos_x_align_plot(n_now(j),m_now(j))*1000 pos_x_align_plot(n_now(j)+dn_mov,m_now(j))*1000],...
                [pos_y_align_plot(n_now(j),m_now(j))*1000 pos_y_align_plot(n_now(j)+dn_mov,m_now(j))*1000],...
                '-','color',cmap_mov(col(n_now(j),m_now(j)),:),'linew',2)
            
            figure(2)
            plot([pos_x_align_plot(n_now(j),m_now(j))*1000 pos_x_align_plot(n_now(j)+dn_mov,m_now(j))*1000],...
                [pos_z_align_plot(n_now(j),m_now(j))*1000 pos_z_align_plot(n_now(j)+dn_mov,m_now(j))*1000],...
                '-','color',cmap_mov(col(n_now(j),m_now(j)),:),'linew',2)
            
            figure(3)
            plot([pos_y_align_plot(n_now(j),m_now(j))*1000 pos_y_align_plot(n_now(j)+dn_mov,m_now(j))*1000],...
                [pos_z_align_plot(n_now(j),m_now(j))*1000 pos_z_align_plot(n_now(j)+dn_mov,m_now(j))*1000],...
                '-','color',cmap_mov(col(n_now(j),m_now(j)),:),'linew',2)
            
            figure(4)
            plot3([pos_x_align_plot(n_now(j),m_now(j))*1000 pos_x_align_plot(n_now(j)+dn_mov,m_now(j))*1000],...
                [pos_y_align_plot(n_now(j),m_now(j))*1000 pos_y_align_plot(n_now(j)+dn_mov,m_now(j))*1000],...
                [pos_z_align_plot(n_now(j),m_now(j))*1000 pos_z_align_plot(n_now(j)+dn_mov,m_now(j))*1000],...
                '-','color',cmap_mov(col(n_now(j),m_now(j)),:),'linew',2)
            
            figure(5)
            plot([pos_x_align_mirror_plot(n_now(j),m_now(j))*1000 pos_x_align_mirror_plot(n_now(j)+dn_mov,m_now(j))*1000],...
                [pos_y_align_plot(n_now(j),m_now(j))*1000 pos_y_align_plot(n_now(j)+dn_mov,m_now(j))*1000],...
                '-','color',cmap_mov(col(n_now(j),m_now(j)),:),'linew',2)
            
            figure(6)
            plot([pos_x_align_mirror_plot(n_now(j),m_now(j))*1000 pos_x_align_mirror_plot(n_now(j)+dn_mov,m_now(j))*1000],...
                [pos_z_align_plot(n_now(j),m_now(j))*1000 pos_z_align_plot(n_now(j)+dn_mov,m_now(j))*1000],...
                '-','color',cmap_mov(col(n_now(j),m_now(j)),:),'linew',2)
            
            figure(7)
            plot([pos_y_align_plot(n_now(j),m_now(j))*1000 pos_y_align_plot(n_now(j)+dn_mov,m_now(j))*1000],...
                [pos_z_align_plot(n_now(j),m_now(j))*1000 pos_z_align_plot(n_now(j)+dn_mov,m_now(j))*1000],...
                '-','color',cmap_mov(col(n_now(j),m_now(j)),:),'linew',2)
            
            figure(8)
            plot3([pos_x_align_mirror_plot(n_now(j),m_now(j))*1000 pos_x_align_mirror_plot(n_now(j)+dn_mov,m_now(j))*1000],...
                [pos_y_align_plot(n_now(j),m_now(j))*1000 pos_y_align_plot(n_now(j)+dn_mov,m_now(j))*1000],...
                [pos_z_align_plot(n_now(j),m_now(j))*1000 pos_z_align_plot(n_now(j)+dn_mov,m_now(j))*1000],...
                '-','color',cmap_mov(col(n_now(j),m_now(j)),:),'linew',2)
        end
    end
    
    if make_mov == 1

        figure(1)
        saveas(gcf,['saccades_topview_TimeColor',num2str(i-mov_start+1),'.tif']);

        figure(2)
        saveas(gcf,['saccades_frontview_TimeColor',num2str(i-mov_start+1),'.tif']);

        figure(3)
        saveas(gcf,['saccades_sideview_TimeColor',num2str(i-mov_start+1),'.tif']);

        figure(4)
        view(-225,30)
        saveas(gcf,['saccades_perspview1_TimeColor',num2str(i-mov_start+1),'.tif']);
        view(-135,30)
        saveas(gcf,['saccades_perspview2_TimeColor',num2str(i-mov_start+1),'.tif']);

        figure(5)
        saveas(gcf,['saccades_topview_TimeColor_mirror',num2str(i-mov_start+1),'.tif']);

        figure(6)
        saveas(gcf,['saccades_frontview_TimeColor_mirror',num2str(i-mov_start+1),'.tif']);

        figure(7)
        saveas(gcf,['saccades_sideview_TimeColor_mirror',num2str(i-mov_start+1),'.tif']);

        figure(8)
        view(-225,30)
        saveas(gcf,['saccades_perspview1_TimeColor_mirror',num2str(i-mov_start+1),'.tif']);
        view(-135,30)
        saveas(gcf,['saccades_perspview2_TimeColor_mirror',num2str(i-mov_start+1),'.tif']);
    end
end

if make_mov == 1
    cd ..
end

figure(1)
axis on
set(gca,'XTick',-10:5:10,'fontsize',20) 
set(gca,'YTick',-10:5:10,'fontsize',20)
xlabel('x [mm]','fontsize',20)
ylabel('y [mm]','fontsize',20)
saveas(gcf,['saccades_topview_TimeColor_time',num2str(t_start),'to',num2str(t_stop),'sec.fig']);
saveas(gcf,['saccades_topview_TimeColor_time',num2str(t_start),'to',num2str(t_stop),'sec.tif']);
plot2svg(['saccades_topview_TimeColor_time',num2str(t_start),'to',num2str(t_stop),'sec.svg']);

figure(2)
axis on
set(gca,'XTick',-10:5:10,'fontsize',20)
set(gca,'YTick',-10:5:10,'fontsize',20)
xlabel('x [mm]','fontsize',20)
ylabel('z [mm]','fontsize',20)
saveas(gcf,['saccades_frontview_TimeColor_time',num2str(t_start),'to',num2str(t_stop),'sec.fig']);
saveas(gcf,['saccades_frontview_TimeColor_time',num2str(t_start),'to',num2str(t_stop),'sec.tif']);
plot2svg(['saccades_frontview_TimeColor_time',num2str(t_start),'to',num2str(t_stop),'sec.svg']);

figure(3)
axis on
set(gca,'XTick',-10:5:10,'fontsize',20) 
set(gca,'YTick',-10:5:10,'fontsize',20)
xlabel('y [mm]','fontsize',20)
ylabel('z [mm]','fontsize',20)
saveas(gcf,['saccades_sideview_TimeColor_time',num2str(t_start),'to',num2str(t_stop),'sec.fig']);
saveas(gcf,['saccades_sideview_TimeColor_time',num2str(t_start),'to',num2str(t_stop),'sec.tif']);
plot2svg(['saccades_sideview_TimeColor_time',num2str(t_start),'to',num2str(t_stop),'sec.svg']);

figure(4)
axis on
set(gca,'XTick',-10:5:10,'fontsize',20) 
set(gca,'YTick',-10:5:10,'fontsize',20)
set(gca,'ZTick',-10:5:10,'fontsize',20)
xlabel('x [mm]','fontsize',20)
ylabel('y [mm]','fontsize',20)
zlabel('z [mm]','fontsize',20)
view(-225,30)
saveas(gcf,['saccades_perspview1_TimeColor_time',num2str(t_start),'to',num2str(t_stop),'sec.fig']);
saveas(gcf,['saccades_perspview1_TimeColor_time',num2str(t_start),'to',num2str(t_stop),'sec.tif']);
plot2svg(['saccades_perspview1_TimeColor_time',num2str(t_start),'to',num2str(t_stop),'sec.svg']);
view(-135,30)
saveas(gcf,['saccades_perspview2_TimeColor_time',num2str(t_start),'to',num2str(t_stop),'sec.fig']);
saveas(gcf,['saccades_perspview2_TimeColor_time',num2str(t_start),'to',num2str(t_stop),'sec.tif']);
plot2svg(['saccades_perspview2_TimeColor_time',num2str(t_start),'to',num2str(t_stop),'sec.svg']);


figure(5)
axis on
set(gca,'XTick',-10:5:10,'fontsize',20) 
set(gca,'YTick',-10:5:10,'fontsize',20)
xlabel('x [mm]','fontsize',20)
ylabel('y [mm]','fontsize',20)
saveas(gcf,['saccades_topview_TimeColor_time',num2str(t_start),'to',num2str(t_stop),'sec_mirror.fig']);
saveas(gcf,['saccades_topview_TimeColor_time',num2str(t_start),'to',num2str(t_stop),'sec_mirror.tif']);
plot2svg(['saccades_topview_TimeColor_time',num2str(t_start),'to',num2str(t_stop),'sec_mirror.svg']);

figure(6)
axis on
set(gca,'XTick',-10:5:10,'fontsize',20)
set(gca,'YTick',-10:5:10,'fontsize',20)
xlabel('x [mm]','fontsize',20)
ylabel('z [mm]','fontsize',20)
saveas(gcf,['saccades_frontview_TimeColor_time',num2str(t_start),'to',num2str(t_stop),'sec_mirror.fig']);
saveas(gcf,['saccades_frontview_TimeColor_time',num2str(t_start),'to',num2str(t_stop),'sec_mirror.tif']);
plot2svg(['saccades_frontview_TimeColor_time',num2str(t_start),'to',num2str(t_stop),'sec_mirror.svg']);

figure(7)
axis on
set(gca,'XTick',-10:5:10,'fontsize',20) 
set(gca,'YTick',-10:5:10,'fontsize',20)
xlabel('y [mm]','fontsize',20)
ylabel('z [mm]','fontsize',20)
saveas(gcf,['saccades_sideview_TimeColor_time',num2str(t_start),'to',num2str(t_stop),'sec_mirror.fig']);
saveas(gcf,['saccades_sideview_TimeColor_time',num2str(t_start),'to',num2str(t_stop),'sec_mirror.tif']);
plot2svg(['saccades_sideview_TimeColor_time',num2str(t_start),'to',num2str(t_stop),'sec_mirror.svg']);

figure(8)
axis on
set(gca,'XTick',-10:5:10,'fontsize',20) 
set(gca,'YTick',-10:5:10,'fontsize',20)
set(gca,'ZTick',-10:5:10,'fontsize',20)
xlabel('x [mm]','fontsize',20)
ylabel('y [mm]','fontsize',20)
zlabel('z [mm]','fontsize',20)
view(-225,30)
saveas(gcf,['saccades_perspview1_TimeColor_time',num2str(t_start),'to',num2str(t_stop),'sec_mirror.fig']);
saveas(gcf,['saccades_perspview1_TimeColor_time',num2str(t_start),'to',num2str(t_stop),'sec_mirror.tif']);
plot2svg(['saccades_perspview1_TimeColor_time',num2str(t_start),'to',num2str(t_stop),'sec_mirror.svg']);
view(-135,30)
saveas(gcf,['saccades_perspview2_TimeColor_time',num2str(t_start),'to',num2str(t_stop),'sec_mirror.fig']);
saveas(gcf,['saccades_perspview2_TimeColor_time',num2str(t_start),'to',num2str(t_stop),'sec_mirror.tif']);
plot2svg(['saccades_perspview2_TimeColor_time',num2str(t_start),'to',num2str(t_stop),'sec_mirror.svg']);

%% movie accelColor
Amax = 10;
% Amax = max(abs(A(:)));

cmap_mov=jet(cmap_steps);
col = round(A./Amax*cmap_steps).*pos_on;
col(col>cmap_steps)=cmap_steps;

figure
colormap(cmap_mov)
colorbar
saveas(gcf,['colorbar_Accel.fig']);
saveas(gcf,['colorbar_Accel.tif']);
plot2svg(['colorbar_Accel.svg']);

close all
figure(1)
hold on
axis equal
axis([-10 10 -10 10])
axis off

figure(2)
hold on
axis equal
axis([-10 10 -5 5])
axis off

figure(3)
hold on
axis equal
axis([-10 10 -5 5])
axis off

figure(4)
hold on
axis equal
axis([-10 10 -10 10 -5 5])
axis off

figure(5)
hold on
axis equal
axis([-1 10 -10 10])
axis off

figure(6)
hold on
axis equal
axis([-1 10 -5 5])
axis off

figure(7)
hold on
axis equal
axis([-10 10 -5 5])
axis off

figure(8)
hold on
axis equal
axis([-1 10 -10 10 -5 5])
axis off

if make_mov == 1
    mkdir('saccades_accelColor')
    cd('saccades_accelColor')

    figure(1)
    saveas(gcf,['saccades_topview_accelColor_Amax',num2str(Amax),'mps2_',num2str(0),'.tif']);

    figure(2)
    saveas(gcf,['saccades_frontview_accelColor_Amax',num2str(Amax),'mps2_',num2str(0),'.tif']);

    figure(3)
    saveas(gcf,['saccades_sideview_accelColor_Amax',num2str(Amax),'mps2_',num2str(0),'.tif']);

    figure(4)
    view(-225,30)
    saveas(gcf,['saccades_perspview1_accelColor_Amax',num2str(Amax),'mps2_',num2str(0),'.tif']);
    view(-135,30)
    saveas(gcf,['saccades_perspview2_accelColor_Amax',num2str(Amax),'mps2_',num2str(0),'.tif']);

    figure(5)
    saveas(gcf,['saccades_topview_accelColor_Amax',num2str(Amax),'mps2_mirror_',num2str(0),'.tif']);

    figure(6)
    saveas(gcf,['saccades_frontview_accelColor_Amax',num2str(Amax),'mps2_mirror_',num2str(0),'.tif']);

    figure(7)
    saveas(gcf,['saccades_sideview_accelColor_Amax',num2str(Amax),'mps2_mirror_',num2str(0),'.tif']);

    figure(8)
    view(-225,30)
    saveas(gcf,['saccades_perspview1_accelColor_Amax',num2str(Amax),'mps2_mirror_',num2str(0),'.tif']);
    view(-135,30)
    saveas(gcf,['saccades_perspview2_accelColor_Amax',num2str(Amax),'mps2_mirror_',num2str(0),'.tif']);
end


if make_mov == 1
    dn_mov = 3;
else
    dn_mov = 40;
end
mov_start = min(n_pos_plot(:));
mov_stop = max(n_pos_plot(:));

for i=mov_start:dn_mov:mov_stop
    
    [n_now,m_now] = find(n_pos_plot==i);

    for j=1:length(n_now)
        if isnan(col(n_now(j),m_now(j)))==0 && isnan(col(n_now(j)+dn_mov,m_now(j)))==0
            
            figure(1)
            plot([pos_x_align_plot(n_now(j),m_now(j))*1000 pos_x_align_plot(n_now(j)+dn_mov,m_now(j))*1000],...
                [pos_y_align_plot(n_now(j),m_now(j))*1000 pos_y_align_plot(n_now(j)+dn_mov,m_now(j))*1000],...
                '-','color',cmap_mov(col(n_now(j),m_now(j)),:),'linew',2)
            
            figure(2)
            plot([pos_x_align_plot(n_now(j),m_now(j))*1000 pos_x_align_plot(n_now(j)+dn_mov,m_now(j))*1000],...
                [pos_z_align_plot(n_now(j),m_now(j))*1000 pos_z_align_plot(n_now(j)+dn_mov,m_now(j))*1000],...
                '-','color',cmap_mov(col(n_now(j),m_now(j)),:),'linew',2)
            
            figure(3)
            plot([pos_y_align_plot(n_now(j),m_now(j))*1000 pos_y_align_plot(n_now(j)+dn_mov,m_now(j))*1000],...
                [pos_z_align_plot(n_now(j),m_now(j))*1000 pos_z_align_plot(n_now(j)+dn_mov,m_now(j))*1000],...
                '-','color',cmap_mov(col(n_now(j),m_now(j)),:),'linew',2)
            
            figure(4)
            plot3([pos_x_align_plot(n_now(j),m_now(j))*1000 pos_x_align_plot(n_now(j)+dn_mov,m_now(j))*1000],...
                [pos_y_align_plot(n_now(j),m_now(j))*1000 pos_y_align_plot(n_now(j)+dn_mov,m_now(j))*1000],...
                [pos_z_align_plot(n_now(j),m_now(j))*1000 pos_z_align_plot(n_now(j)+dn_mov,m_now(j))*1000],...
                '-','color',cmap_mov(col(n_now(j),m_now(j)),:),'linew',2)
            
            figure(5)
            plot([pos_x_align_mirror_plot(n_now(j),m_now(j))*1000 pos_x_align_mirror_plot(n_now(j)+dn_mov,m_now(j))*1000],...
                [pos_y_align_plot(n_now(j),m_now(j))*1000 pos_y_align_plot(n_now(j)+dn_mov,m_now(j))*1000],...
                '-','color',cmap_mov(col(n_now(j),m_now(j)),:),'linew',2)
            
            figure(6)
            plot([pos_x_align_mirror_plot(n_now(j),m_now(j))*1000 pos_x_align_mirror_plot(n_now(j)+dn_mov,m_now(j))*1000],...
                [pos_z_align_plot(n_now(j),m_now(j))*1000 pos_z_align_plot(n_now(j)+dn_mov,m_now(j))*1000],...
                '-','color',cmap_mov(col(n_now(j),m_now(j)),:),'linew',2)
            
            figure(7)
            plot([pos_y_align_plot(n_now(j),m_now(j))*1000 pos_y_align_plot(n_now(j)+dn_mov,m_now(j))*1000],...
                [pos_z_align_plot(n_now(j),m_now(j))*1000 pos_z_align_plot(n_now(j)+dn_mov,m_now(j))*1000],...
                '-','color',cmap_mov(col(n_now(j),m_now(j)),:),'linew',2)
            
            figure(8)
            plot3([pos_x_align_mirror_plot(n_now(j),m_now(j))*1000 pos_x_align_mirror_plot(n_now(j)+dn_mov,m_now(j))*1000],...
                [pos_y_align_plot(n_now(j),m_now(j))*1000 pos_y_align_plot(n_now(j)+dn_mov,m_now(j))*1000],...
                [pos_z_align_plot(n_now(j),m_now(j))*1000 pos_z_align_plot(n_now(j)+dn_mov,m_now(j))*1000],...
                '-','color',cmap_mov(col(n_now(j),m_now(j)),:),'linew',2)
        end
    end
    
    if make_mov == 1

        figure(1)
        saveas(gcf,['saccades_topview_accelColor_Amax',num2str(Amax),'mps2_',num2str(i-mov_start+1),'.tif']);

        figure(2)
        saveas(gcf,['saccades_frontview_accelColor_Amax',num2str(Amax),'mps2_',num2str(i-mov_start+1),'.tif']);

        figure(3)
        saveas(gcf,['saccades_sideview_accelColor_Amax',num2str(Amax),'mps2_',num2str(i-mov_start+1),'.tif']);

        figure(4)
        view(-225,30)
        saveas(gcf,['saccades_perspview1_accelColor_Amax',num2str(Amax),'mps2_',num2str(i-mov_start+1),'.tif']);
        view(-135,30)
        saveas(gcf,['saccades_perspview2_accelColor_Amax',num2str(Amax),'mps2_',num2str(i-mov_start+1),'.tif']);

        figure(5)
        saveas(gcf,['saccades_topview_accelColor_Amax',num2str(Amax),'mps2_mirror_',num2str(i-mov_start+1),'.tif']);

        figure(6)
        saveas(gcf,['saccades_frontview_accelColor_Amax',num2str(Amax),'mps2_mirror_',num2str(i-mov_start+1),'.tif']);

        figure(7)
        saveas(gcf,['saccades_sideview_accelColor_Amax',num2str(Amax),'mps2_mirror_',num2str(i-mov_start+1),'.tif']);

        figure(8)
        view(-225,30)
        saveas(gcf,['saccades_perspview1_accelColor_Amax',num2str(Amax),'mps2_mirror_',num2str(i-mov_start+1),'.tif']);
        view(-135,30)
        saveas(gcf,['saccades_perspview2_accelColor_Amax',num2str(Amax),'mps2_mirror_',num2str(i-mov_start+1),'.tif']);
    end
end

if make_mov == 1
    cd ..
end

figure(1)
axis on
set(gca,'XTick',-10:5:10,'fontsize',20) 
set(gca,'YTick',-10:5:10,'fontsize',20)
xlabel('x [mm]','fontsize',20)
ylabel('y [mm]','fontsize',20)
saveas(gcf,['saccades_topview_accelColor_Amax',num2str(Amax),'mps2_time',num2str(t_start),'to',num2str(t_stop),'sec.fig']);
saveas(gcf,['saccades_topview_accelColor_Amax',num2str(Amax),'mps2_time',num2str(t_start),'to',num2str(t_stop),'sec.tif']);
plot2svg(['saccades_topview_accelColor_Amax',num2str(Amax),'mps2_time',num2str(t_start),'to',num2str(t_stop),'sec.svg']);

figure(2)
axis on
set(gca,'XTick',-10:5:10,'fontsize',20) 
set(gca,'YTick',-10:5:10,'fontsize',20)
xlabel('x [mm]','fontsize',20)
ylabel('z [mm]','fontsize',20)
saveas(gcf,['saccades_frontview_accelColor_Amax',num2str(Amax),'mps2_time',num2str(t_start),'to',num2str(t_stop),'sec.fig']);
saveas(gcf,['saccades_frontview_accelColor_Amax',num2str(Amax),'mps2_time',num2str(t_start),'to',num2str(t_stop),'sec.tif']);
plot2svg(['saccades_frontview_accelColor_Amax',num2str(Amax),'mps2_time',num2str(t_start),'to',num2str(t_stop),'sec.svg']);

figure(3)
axis on
set(gca,'XTick',-10:5:10,'fontsize',20) 
set(gca,'YTick',-10:5:10,'fontsize',20)
xlabel('y [mm]','fontsize',20)
ylabel('z [mm]','fontsize',20)
saveas(gcf,['saccades_sideview_accelColor_Amax',num2str(Amax),'mps2_time',num2str(t_start),'to',num2str(t_stop),'sec.fig']);
saveas(gcf,['saccades_sideview_accelColor_Amax',num2str(Amax),'mps2_time',num2str(t_start),'to',num2str(t_stop),'sec.tif']);
plot2svg(['saccades_sideview_accelColor_Amax',num2str(Amax),'mps2_time',num2str(t_start),'to',num2str(t_stop),'sec.svg']);

figure(4)
axis on
set(gca,'XTick',-10:5:10,'fontsize',20) 
set(gca,'YTick',-10:5:10,'fontsize',20)
set(gca,'ZTick',-10:5:10,'fontsize',20)
xlabel('x [mm]','fontsize',20)
ylabel('y [mm]','fontsize',20)
zlabel('z [mm]','fontsize',20)
view(-225,30)
saveas(gcf,['saccades_perspview1_accelColor_Amax',num2str(Amax),'mps2_time',num2str(t_start),'to',num2str(t_stop),'sec.fig']);
saveas(gcf,['saccades_perspview1_accelColor_Amax',num2str(Amax),'mps2_time',num2str(t_start),'to',num2str(t_stop),'sec.tif']);
plot2svg(['saccades_perspview1_accelColor_Amax',num2str(Amax),'mps2_time',num2str(t_start),'to',num2str(t_stop),'sec.svg']);
view(-135,30)
saveas(gcf,['saccades_perspview2_accelColor_Amax',num2str(Amax),'mps2_time',num2str(t_start),'to',num2str(t_stop),'sec.fig']);
saveas(gcf,['saccades_perspview2_accelColor_Amax',num2str(Amax),'mps2_time',num2str(t_start),'to',num2str(t_stop),'sec.tif']);
plot2svg(['saccades_perspview2_accelColor_Amax',num2str(Amax),'mps2_time',num2str(t_start),'to',num2str(t_stop),'sec.svg']);


figure(5)
axis on
set(gca,'XTick',-10:5:10,'fontsize',20) 
set(gca,'YTick',-10:5:10,'fontsize',20)
xlabel('x [mm]','fontsize',20)
ylabel('y [mm]','fontsize',20)
saveas(gcf,['saccades_topview_accelColor_Amax',num2str(Amax),'mps2_time',num2str(t_start),'to',num2str(t_stop),'sec_mirror.fig']);
saveas(gcf,['saccades_topview_accelColor_Amax',num2str(Amax),'mps2_time',num2str(t_start),'to',num2str(t_stop),'sec_mirror.tif']);
plot2svg(['saccades_topview_accelColor_Amax',num2str(Amax),'mps2_time',num2str(t_start),'to',num2str(t_stop),'sec_mirror.svg']);

figure(6)
axis on
set(gca,'XTick',-10:5:10,'fontsize',20) 
set(gca,'YTick',-10:5:10,'fontsize',20)
xlabel('x [mm]','fontsize',20)
ylabel('z [mm]','fontsize',20)
saveas(gcf,['saccades_frontview_accelColor_Amax',num2str(Amax),'mps2_time',num2str(t_start),'to',num2str(t_stop),'sec_mirror.fig']);
saveas(gcf,['saccades_frontview_accelColor_Amax',num2str(Amax),'mps2_time',num2str(t_start),'to',num2str(t_stop),'sec_mirror.tif']);
plot2svg(['saccades_frontview_accelColor_Amax',num2str(Amax),'mps2_time',num2str(t_start),'to',num2str(t_stop),'sec_mirror.svg']);

figure(7)
axis on
set(gca,'XTick',-10:5:10,'fontsize',20) 
set(gca,'YTick',-10:5:10,'fontsize',20)
xlabel('y [mm]','fontsize',20)
ylabel('z [mm]','fontsize',20)
saveas(gcf,['saccades_sideview_accelColor_Amax',num2str(Amax),'mps2_time',num2str(t_start),'to',num2str(t_stop),'sec_mirror.fig']);
saveas(gcf,['saccades_sideview_accelColor_Amax',num2str(Amax),'mps2_time',num2str(t_start),'to',num2str(t_stop),'sec_mirror.tif']);
plot2svg(['saccades_sideview_accelColor_Amax',num2str(Amax),'mps2_time',num2str(t_start),'to',num2str(t_stop),'sec_mirror.svg']);

figure(8)
axis on
set(gca,'XTick',-10:5:10,'fontsize',20) 
set(gca,'YTick',-10:5:10,'fontsize',20)
set(gca,'ZTick',-10:5:10,'fontsize',20)
xlabel('x [mm]','fontsize',20)
ylabel('y [mm]','fontsize',20)
zlabel('z [mm]','fontsize',20)
view(-225,30)
saveas(gcf,['saccades_perspview1_accelColor_Amax',num2str(Amax),'mps2_time',num2str(t_start),'to',num2str(t_stop),'sec_mirror.fig']);
saveas(gcf,['saccades_perspview1_accelColor_Amax',num2str(Amax),'mps2_time',num2str(t_start),'to',num2str(t_stop),'sec_mirror.tif']);
plot2svg(['saccades_perspview1_accelColor_Amax',num2str(Amax),'mps2_time',num2str(t_start),'to',num2str(t_stop),'sec_mirror.svg']);
view(-135,30)
saveas(gcf,['saccades_perspview2_accelColor_Amax',num2str(Amax),'mps2_time',num2str(t_start),'to',num2str(t_stop),'sec_mirror.fig']);
saveas(gcf,['saccades_perspview2_accelColor_Amax',num2str(Amax),'mps2_time',num2str(t_start),'to',num2str(t_stop),'sec_mirror.tif']);
plot2svg(['saccades_perspview2_accelColor_Amax',num2str(Amax),'mps2_time',num2str(t_start),'to',num2str(t_stop),'sec_mirror.svg']);

%% movie VelColor
Vmax = 0.3;

cmap_mov=jet(cmap_steps);
col = round(V./Vmax*cmap_steps).*pos_on;
col(col>cmap_steps)=cmap_steps;
    
figure
colormap(cmap_mov)
colorbar
saveas(gcf,['colorbar_Vel.tif']);
saveas(gcf,['colorbar_Vel.fig']);
plot2svg(['colorbar_Vel.svg']);

close all
figure(1)
hold on
axis equal
axis([-10 10 -10 10])
axis off

figure(2)
hold on
axis equal
axis([-10 10 -5 5])
axis off

figure(3)
hold on
axis equal
axis([-10 10 -5 5])
axis off

figure(4)
hold on
axis equal
axis([-10 10 -10 10 -5 5])
axis off

figure(5)
hold on
axis equal
axis([-1 10 -10 10])
axis off

figure(6)
hold on
axis equal
axis([-1 10 -5 5])
axis off

figure(7)
hold on
axis equal
axis([-10 10 -5 5])
axis off

figure(8)
hold on
axis equal
axis([-1 10 -10 10 -5 5])
axis off

if make_mov == 1
    mkdir('saccades_VelColor')
    cd('saccades_VelColor')

    figure(1)
    saveas(gcf,['saccades_topview_VelColor_Amax',num2str(Amax),'mps2_',num2str(0),'.tif']);

    figure(2)
    saveas(gcf,['saccades_frontview_VelColor_Amax',num2str(Amax),'mps2_',num2str(0),'.tif']);

    figure(3)
    saveas(gcf,['saccades_sideview_VelColor_Amax',num2str(Amax),'mps2_',num2str(0),'.tif']);

    figure(4)
    view(-225,30)
    saveas(gcf,['saccades_perspview1_VelColor_Amax',num2str(Amax),'mps2_',num2str(0),'.tif']);
    view(-135,30)
    saveas(gcf,['saccades_perspview2_VelColor_Amax',num2str(Amax),'mps2_',num2str(0),'.tif']);

    figure(5)
    saveas(gcf,['saccades_topview_VelColor_Amax',num2str(Amax),'mps2_mirror_',num2str(0),'.tif']);

    figure(6)
    saveas(gcf,['saccades_frontview_VelColor_Amax',num2str(Amax),'mps2_mirror_',num2str(0),'.tif']);

    figure(7)
    saveas(gcf,['saccades_sideview_VelColor_Amax',num2str(Amax),'mps2_mirror_',num2str(0),'.tif']);

    figure(8)
    view(-225,30)
    saveas(gcf,['saccades_perspview1_VelColor_Amax',num2str(Amax),'mps2_mirror_',num2str(0),'.tif']);
    view(-135,30)
    saveas(gcf,['saccades_perspview2_VelColor_Amax',num2str(Amax),'mps2_mirror_',num2str(0),'.tif']);
end

if make_mov == 1
    dn_mov = 3;
else
    dn_mov = 40;
end
mov_start = min(n_pos_plot(:));
mov_stop = max(n_pos_plot(:));

for i=mov_start:dn_mov:mov_stop
    
    [n_now,m_now] = find(n_pos_plot==i);

    for j=1:length(n_now)
        if isnan(col(n_now(j),m_now(j)))==0 && isnan(col(n_now(j)+dn_mov,m_now(j)))==0
            
            figure(1)
            plot([pos_x_align_plot(n_now(j),m_now(j))*1000 pos_x_align_plot(n_now(j)+dn_mov,m_now(j))*1000],...
                [pos_y_align_plot(n_now(j),m_now(j))*1000 pos_y_align_plot(n_now(j)+dn_mov,m_now(j))*1000],...
                '-','color',cmap_mov(col(n_now(j),m_now(j)),:),'linew',2)
            
            figure(2)
            plot([pos_x_align_plot(n_now(j),m_now(j))*1000 pos_x_align_plot(n_now(j)+dn_mov,m_now(j))*1000],...
                [pos_z_align_plot(n_now(j),m_now(j))*1000 pos_z_align_plot(n_now(j)+dn_mov,m_now(j))*1000],...
                '-','color',cmap_mov(col(n_now(j),m_now(j)),:),'linew',2)
            
            figure(3)
            plot([pos_y_align_plot(n_now(j),m_now(j))*1000 pos_y_align_plot(n_now(j)+dn_mov,m_now(j))*1000],...
                [pos_z_align_plot(n_now(j),m_now(j))*1000 pos_z_align_plot(n_now(j)+dn_mov,m_now(j))*1000],...
                '-','color',cmap_mov(col(n_now(j),m_now(j)),:),'linew',2)
            
            figure(4)
            plot3([pos_x_align_plot(n_now(j),m_now(j))*1000 pos_x_align_plot(n_now(j)+dn_mov,m_now(j))*1000],...
                [pos_y_align_plot(n_now(j),m_now(j))*1000 pos_y_align_plot(n_now(j)+dn_mov,m_now(j))*1000],...
                [pos_z_align_plot(n_now(j),m_now(j))*1000 pos_z_align_plot(n_now(j)+dn_mov,m_now(j))*1000],...
                '-','color',cmap_mov(col(n_now(j),m_now(j)),:),'linew',2)
            
            figure(5)
            plot([pos_x_align_mirror_plot(n_now(j),m_now(j))*1000 pos_x_align_mirror_plot(n_now(j)+dn_mov,m_now(j))*1000],...
                [pos_y_align_plot(n_now(j),m_now(j))*1000 pos_y_align_plot(n_now(j)+dn_mov,m_now(j))*1000],...
                '-','color',cmap_mov(col(n_now(j),m_now(j)),:),'linew',2)
            
            figure(6)
            plot([pos_x_align_mirror_plot(n_now(j),m_now(j))*1000 pos_x_align_mirror_plot(n_now(j)+dn_mov,m_now(j))*1000],...
                [pos_z_align_plot(n_now(j),m_now(j))*1000 pos_z_align_plot(n_now(j)+dn_mov,m_now(j))*1000],...
                '-','color',cmap_mov(col(n_now(j),m_now(j)),:),'linew',2)
            
            figure(7)
            plot([pos_y_align_plot(n_now(j),m_now(j))*1000 pos_y_align_plot(n_now(j)+dn_mov,m_now(j))*1000],...
                [pos_z_align_plot(n_now(j),m_now(j))*1000 pos_z_align_plot(n_now(j)+dn_mov,m_now(j))*1000],...
                '-','color',cmap_mov(col(n_now(j),m_now(j)),:),'linew',2)
            
            figure(8)
            plot3([pos_x_align_mirror_plot(n_now(j),m_now(j))*1000 pos_x_align_mirror_plot(n_now(j)+dn_mov,m_now(j))*1000],...
                [pos_y_align_plot(n_now(j),m_now(j))*1000 pos_y_align_plot(n_now(j)+dn_mov,m_now(j))*1000],...
                [pos_z_align_plot(n_now(j),m_now(j))*1000 pos_z_align_plot(n_now(j)+dn_mov,m_now(j))*1000],...
                '-','color',cmap_mov(col(n_now(j),m_now(j)),:),'linew',2)
        end
    end
    
    if make_mov == 1

        figure(1)
        saveas(gcf,['saccades_topview_VelColor_Amax',num2str(Amax),'mps2_',num2str(i-mov_start+1),'.tif']);

        figure(2)
        saveas(gcf,['saccades_frontview_VelColor_Amax',num2str(Amax),'mps2_',num2str(i-mov_start+1),'.tif']);

        figure(3)
        saveas(gcf,['saccades_sideview_VelColor_Amax',num2str(Amax),'mps2_',num2str(i-mov_start+1),'.tif']);

        figure(4)
        view(-225,30)
        saveas(gcf,['saccades_perspview1_VelColor_Amax',num2str(Amax),'mps2_',num2str(i-mov_start+1),'.tif']);
        view(-135,30)
        saveas(gcf,['saccades_perspview2_VelColor_Amax',num2str(Amax),'mps2_',num2str(i-mov_start+1),'.tif']);

        figure(5)
        saveas(gcf,['saccades_topview_VelColor_Amax',num2str(Amax),'mps2_mirror_',num2str(i-mov_start+1),'.tif']);

        figure(6)
        saveas(gcf,['saccades_frontview_VelColor_Amax',num2str(Amax),'mps2_mirror_',num2str(i-mov_start+1),'.tif']);

        figure(7)
        saveas(gcf,['saccades_sideview_VelColor_Amax',num2str(Amax),'mps2_mirror_',num2str(i-mov_start+1),'.tif']);

        figure(8)
        view(-225,30)
        saveas(gcf,['saccades_perspview1_VelColor_Amax',num2str(Amax),'mps2_mirror_',num2str(i-mov_start+1),'.tif']);
        view(-135,30)
        saveas(gcf,['saccades_perspview2_VelColor_Amax',num2str(Amax),'mps2_mirror_',num2str(i-mov_start+1),'.tif']);
    end
end

if make_mov == 1
    cd ..
end

figure(1)
axis on
set(gca,'XTick',-10:5:10,'fontsize',20) 
set(gca,'YTick',-10:5:10,'fontsize',20)
xlabel('x [mm]','fontsize',20)
ylabel('y [mm]','fontsize',20)
saveas(gcf,['saccades_topview_VelColor_Amax',num2str(Amax),'mps2_time',num2str(t_start),'to',num2str(t_stop),'sec.fig']);
saveas(gcf,['saccades_topview_VelColor_Amax',num2str(Amax),'mps2_time',num2str(t_start),'to',num2str(t_stop),'sec.tif']);
plot2svg(['saccades_topview_VelColor_Amax',num2str(Amax),'mps2_time',num2str(t_start),'to',num2str(t_stop),'sec.svg']);

figure(2)
axis on
set(gca,'XTick',-10:5:10,'fontsize',20) 
set(gca,'YTick',-10:5:10,'fontsize',20)
xlabel('x [mm]','fontsize',20)
ylabel('z [mm]','fontsize',20)
saveas(gcf,['saccades_frontview_VelColor_Amax',num2str(Amax),'mps2_time',num2str(t_start),'to',num2str(t_stop),'sec.fig']);
saveas(gcf,['saccades_frontview_VelColor_Amax',num2str(Amax),'mps2_time',num2str(t_start),'to',num2str(t_stop),'sec.tif']);
plot2svg(['saccades_frontview_VelColor_Amax',num2str(Amax),'mps2_time',num2str(t_start),'to',num2str(t_stop),'sec.svg']);

figure(3)
axis on
set(gca,'XTick',-10:5:10,'fontsize',20) 
set(gca,'YTick',-10:5:10,'fontsize',20)
xlabel('y [mm]','fontsize',20)
ylabel('z [mm]','fontsize',20)
saveas(gcf,['saccades_sideview_VelColor_Amax',num2str(Amax),'mps2_time',num2str(t_start),'to',num2str(t_stop),'sec.fig']);
saveas(gcf,['saccades_sideview_VelColor_Amax',num2str(Amax),'mps2_time',num2str(t_start),'to',num2str(t_stop),'sec.tif']);
plot2svg(['saccades_sideview_VelColor_Amax',num2str(Amax),'mps2_time',num2str(t_start),'to',num2str(t_stop),'sec.svg']);

figure(4)
axis on
set(gca,'XTick',-10:5:10,'fontsize',20) 
set(gca,'YTick',-10:5:10,'fontsize',20)
set(gca,'ZTick',-10:5:10,'fontsize',20)
xlabel('x [mm]','fontsize',20)
ylabel('y [mm]','fontsize',20)
zlabel('z [mm]','fontsize',20)
view(-225,30)
saveas(gcf,['saccades_perspview1_VelColor_Amax',num2str(Amax),'mps2_time',num2str(t_start),'to',num2str(t_stop),'sec.fig']);
saveas(gcf,['saccades_perspview1_VelColor_Amax',num2str(Amax),'mps2_time',num2str(t_start),'to',num2str(t_stop),'sec.tif']);
plot2svg(['saccades_perspview1_VelColor_Amax',num2str(Amax),'mps2_time',num2str(t_start),'to',num2str(t_stop),'sec.svg']);
view(-135,30)
saveas(gcf,['saccades_perspview2_VelColor_Amax',num2str(Amax),'mps2_time',num2str(t_start),'to',num2str(t_stop),'sec.fig']);
saveas(gcf,['saccades_perspview2_VelColor_Amax',num2str(Amax),'mps2_time',num2str(t_start),'to',num2str(t_stop),'sec.tif']);
plot2svg(['saccades_perspview2_VelColor_Amax',num2str(Amax),'mps2_time',num2str(t_start),'to',num2str(t_stop),'sec.svg']);


figure(5)
axis on
set(gca,'XTick',-10:5:10,'fontsize',20) 
set(gca,'YTick',-10:5:10,'fontsize',20)
xlabel('x [mm]','fontsize',20)
ylabel('y [mm]','fontsize',20)
saveas(gcf,['saccades_topview_VelColor_Amax',num2str(Amax),'mps2_time',num2str(t_start),'to',num2str(t_stop),'sec_mirror.fig']);
saveas(gcf,['saccades_topview_VelColor_Amax',num2str(Amax),'mps2_time',num2str(t_start),'to',num2str(t_stop),'sec_mirror.tif']);
plot2svg(['saccades_topview_VelColor_Amax',num2str(Amax),'mps2_time',num2str(t_start),'to',num2str(t_stop),'sec_mirror.svg']);

figure(6)
axis on
set(gca,'XTick',-10:5:10,'fontsize',20) 
set(gca,'YTick',-10:5:10,'fontsize',20)
xlabel('x [mm]','fontsize',20)
ylabel('z [mm]','fontsize',20)
saveas(gcf,['saccades_frontview_VelColor_Amax',num2str(Amax),'mps2_time',num2str(t_start),'to',num2str(t_stop),'sec_mirror.fig']);
saveas(gcf,['saccades_frontview_VelColor_Amax',num2str(Amax),'mps2_time',num2str(t_start),'to',num2str(t_stop),'sec_mirror.tif']);
plot2svg(['saccades_frontview_VelColor_Amax',num2str(Amax),'mps2_time',num2str(t_start),'to',num2str(t_stop),'sec_mirror.svg']);

figure(7)
axis on
set(gca,'XTick',-10:5:10,'fontsize',20) 
set(gca,'YTick',-10:5:10,'fontsize',20)
xlabel('y [mm]','fontsize',20)
ylabel('z [mm]','fontsize',20)
saveas(gcf,['saccades_sideview_VelColor_Amax',num2str(Amax),'mps2_time',num2str(t_start),'to',num2str(t_stop),'sec_mirror.fig']);
saveas(gcf,['saccades_sideview_VelColor_Amax',num2str(Amax),'mps2_time',num2str(t_start),'to',num2str(t_stop),'sec_mirror.tif']);
plot2svg(['saccades_sideview_VelColor_Amax',num2str(Amax),'mps2_time',num2str(t_start),'to',num2str(t_stop),'sec_mirror.svg']);

figure(8)
axis on
set(gca,'XTick',-10:5:10,'fontsize',20) 
set(gca,'YTick',-10:5:10,'fontsize',20)
set(gca,'ZTick',-10:5:10,'fontsize',20)
xlabel('x [mm]','fontsize',20)
ylabel('y [mm]','fontsize',20)
zlabel('z [mm]','fontsize',20)
view(-225,30)
saveas(gcf,['saccades_perspview1_VelColor_Amax',num2str(Amax),'mps2_time',num2str(t_start),'to',num2str(t_stop),'sec_mirror.fig']);
saveas(gcf,['saccades_perspview1_VelColor_Amax',num2str(Amax),'mps2_time',num2str(t_start),'to',num2str(t_stop),'sec_mirror.tif']);
plot2svg(['saccades_perspview1_VelColor_Amax',num2str(Amax),'mps2_time',num2str(t_start),'to',num2str(t_stop),'sec_mirror.svg']);
view(-135,30)
saveas(gcf,['saccades_perspview2_VelColor_Amax',num2str(Amax),'mps2_time',num2str(t_start),'to',num2str(t_stop),'sec_mirror.fig']);
saveas(gcf,['saccades_perspview2_VelColor_Amax',num2str(Amax),'mps2_time',num2str(t_start),'to',num2str(t_stop),'sec_mirror.tif']);
plot2svg(['saccades_perspview2_VelColor_Amax',num2str(Amax),'mps2_time',num2str(t_start),'to',num2str(t_stop),'sec_mirror.svg']);

cd ..

end