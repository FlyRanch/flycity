function body_and_wing_model_plots( settings, pathDB )

seq_nr = 98;

R_strk   = pathDB.rot_mat.Rstr;

body_pos = [0; 0; 0];
Rb       = R_strk'*[1 0 0; 0 -1 0; 0 0 -1];
RL       = [1 0 0; 0 cos(pi/3) sin(pi/3) ; 0 -sin(pi/3) cos(pi/3)];
RR       = [1 0 0; 0 cos(-pi/3) sin(-pi/3) ; 0 -sin(-pi/3) cos(-pi/3)];

body_model.x_mod        = pathDB.body_model.x_mod(:,:,seq_nr);
body_model.y_mod        = pathDB.body_model.y_mod(:,:,seq_nr);
body_model.z_mod        = pathDB.body_model.z_mod(:,:,seq_nr);
body_model.R_strk       = R_strk;

wing_model.length       = pathDB.wing_model.length(seq_nr);

wing_model.x_mod_L      = pathDB.wing_model.x_mod_L(:,:,seq_nr);
wing_model.y_mod_L      = pathDB.wing_model.y_mod_L(:,:,seq_nr);
wing_model.z_mod_L      = pathDB.wing_model.z_mod_L(:,:,seq_nr);
    
wing_model.x_mod_R      = pathDB.wing_model.x_mod_R(:,:,seq_nr);
wing_model.y_mod_R      = pathDB.wing_model.y_mod_R(:,:,seq_nr);
wing_model.z_mod_R      = pathDB.wing_model.z_mod_R(:,:,seq_nr);
    
body_model.Joint_left   = pathDB.body_model.Joint_left(seq_nr,:)';
body_model.Joint_right  = pathDB.body_model.Joint_right(seq_nr,:)';

figure(1)
hFig = figure(1);
set(gcf,'PaperPositionMode','auto');
set(hFig,'Position',[0 0 800 800]);
hold on
Fly_model_plot( body_pos, Rb, RL, RR, body_model, wing_model )
title('Body and wing reference frames')
xlabel('x [mm]')
ylabel('y [mm]')
zlabel('z [mm]')
xlim([-3 3])
ylim([-3 3])
zlim([-3 3])
hold off
saveas(hFig,[char(settings.plot_loc) '/body_wing_model/axes_model'],'fig')

sect_nr_w_x = 8;
sect_nr_w_y = 20;

J_pos_R = body_model.Joint_right;

x_cont_w = -J_pos_R(1) + wing_model.x_mod_R(:,1);
y_cont_w = -J_pos_R(2) + wing_model.y_mod_R(:,1);

clear temp_x_cont_b temp_y_cont_b temp_x_cont_w temp_y_cont_w

w_length = wing_model.length;
w_width = max(x_cont_w)-min(x_cont_w);

delta_x_w = w_width/sect_nr_w_x;
delta_y_w = w_length/sect_nr_w_y;

pos_w_x = min(x_cont_w)+((0.5/(sect_nr_w_x)):(1/(sect_nr_w_x)):(1-(0.5/(sect_nr_w_x)))).*w_width;
pos_w_y = ((0.5/(sect_nr_w_y)):(1/(sect_nr_w_y)):(1-(0.5/(sect_nr_w_y)))).*w_length;
pos_w_z = ((-0.5+(0.5/(sect_nr_w_x))):(1/(sect_nr_w_x)):(1-(0.5/(sect_nr_w_x)))).*w_width;

[Xw, Yw] = meshgrid(pos_w_x, pos_w_y);

figure(2)
hFig = figure(2);
set(gcf,'PaperPositionMode','auto');
set(hFig,'Position',[0 0 1200 800]);
hold on
surf(J_pos_R(2)+wing_model.y_mod_L,-J_pos_R(1)+wing_model.x_mod_L,-J_pos_R(3)+wing_model.z_mod_L,'facecolor',[0.3 0.3 0.3],'edgecolor','none')
plot(-Yw,Xw,'o','MarkerSize',1.5,'MarkerEdgeColor', 'b')
xlim([-3.1 0])
ylim([-0.9 0.5])
title('Wing raster 2D')
xlabel('y [mm]')
ylabel('x [mm]')
axis equal
hold off
saveas(hFig,[char(settings.plot_loc) '/body_wing_model/wing_2D_raster'],'fig')

theta_L2    = (10/180)*pi;
eta_L2      = -(135/180)*pi;
phi_L2      = -(45/180)*pi;

theta_R2    = -(10/180)*pi;
eta_R2      = -(135/180)*pi;
phi_R2      = (45/180)*pi;

RL2 =   [-1 0 0; 0 1 0; 0 0 -1]*[cos(eta_L2) 0 -sin(eta_L2); 0 1 0; sin(eta_L2) 0 cos(eta_L2)]*...
        [1 0 0; 0 cos(theta_L2) sin(theta_L2); 0 -sin(theta_L2) cos(theta_L2)]*...
        [cos(phi_L2) sin(phi_L2) 0; -sin(phi_L2) cos(phi_L2) 0; 0 0 1]*R_strk;
RR2 =   [-1 0 0; 0 1 0; 0 0 -1]*[cos(eta_R2) 0 -sin(eta_R2); 0 1 0; sin(eta_R2) 0 cos(eta_R2)]*...
        [1 0 0; 0 cos(theta_R2) sin(theta_R2); 0 -sin(theta_R2) cos(theta_R2)]*...
        [cos(phi_R2) sin(phi_R2) 0; -sin(phi_R2) cos(phi_R2) 0; 0 0 1]*R_strk;
    
theta_L3    = (25/180)*pi;
eta_L3      = -(135/180)*pi;
phi_L3      = -(20/180)*pi;

theta_R3    = -(25/180)*pi;
eta_R3      = -(135/180)*pi;
phi_R3      = (20/180)*pi;

RL3 =   [-1 0 0; 0 1 0; 0 0 -1]*[cos(eta_L3) 0 -sin(eta_L3); 0 1 0; sin(eta_L3) 0 cos(eta_L3)]*...
        [1 0 0; 0 cos(theta_L3) sin(theta_L3); 0 -sin(theta_L3) cos(theta_L3)]*...
        [cos(phi_L3) sin(phi_L3) 0; -sin(phi_L3) cos(phi_L3) 0; 0 0 1]*R_strk;
RR3 =   [-1 0 0; 0 1 0; 0 0 -1]*[cos(eta_R3) 0 -sin(eta_R3); 0 1 0; sin(eta_R3) 0 cos(eta_R3)]*...
        [1 0 0; 0 cos(theta_R3) sin(theta_R3); 0 -sin(theta_R3) cos(theta_R3)]*...
        [cos(phi_R2) sin(phi_R3) 0; -sin(phi_R3) cos(phi_R3) 0; 0 0 1]*R_strk;
    
theta_L4    = (0/180)*pi;
eta_L4      = -(45/180)*pi;
phi_L4      = -(20/180)*pi;

theta_R4    = -(0/180)*pi;
eta_R4      = -(45/180)*pi;
phi_R4      = (20/180)*pi;

RL4 =   [-1 0 0; 0 1 0; 0 0 -1]*[cos(eta_L4) 0 -sin(eta_L4); 0 1 0; sin(eta_L4) 0 cos(eta_L4)]*...
        [1 0 0; 0 cos(theta_L4) sin(theta_L4); 0 -sin(theta_L4) cos(theta_L4)]*...
        [cos(phi_L4) sin(phi_L4) 0; -sin(phi_L4) cos(phi_L4) 0; 0 0 1]*R_strk;
RR4 =   [-1 0 0; 0 1 0; 0 0 -1]*[cos(eta_R4) 0 -sin(eta_R4); 0 1 0; sin(eta_R4) 0 cos(eta_R4)]*...
        [1 0 0; 0 cos(theta_R4) sin(theta_R4); 0 -sin(theta_R4) cos(theta_R4)]*...
        [cos(phi_R4) sin(phi_R4) 0; -sin(phi_R4) cos(phi_R4) 0; 0 0 1]*R_strk;

figure(3)
hFig = figure(3);
set(gcf,'PaperPositionMode','auto');
set(hFig,'Position',[0 0 1200 800]);
hold on
subplot(2,2,[1; 3]); strokeplane_ref_frame_plot( Rb, RL2, RR2, body_model, wing_model )
xlabel('x [mm]')
ylabel('y [mm]')
zlabel('z [mm]')
subplot(2,2,2); strokeplane_ref_frame_plot( Rb, RL2, RR2, body_model, wing_model )
xlabel('x [mm]')
ylabel('y [mm]')
zlabel('z [mm]')
subplot(2,2,4); strokeplane_ref_frame_plot( Rb, RL4, RR4, body_model, wing_model )
xlabel('x [mm]')
ylabel('y [mm]')
zlabel('z [mm]')
hold off
[~,h1] = suplabel('Strokeplane Reference Frame', 't');
set(h1,'FontSize',10)
saveas(hFig,[char(settings.plot_loc) '/body_wing_model/strokeplane_ref_frame'],'fig')

figure(4)
hFig = figure(4);
set(gcf,'PaperPositionMode','auto');
set(hFig,'Position',[0 0 800 800]);
hold on
Blade_element_plot( Rb, RL2, RR2, body_model, wing_model )
title('Blade element lay-out')
xlabel('x [mm]')
ylabel('y [mm]')
zlabel('z [mm]')
hold off
saveas(hFig,[char(settings.plot_loc) '/body_wing_model/blade_element'],'fig')

end

% % [Xw2,Yw2,Zw] = meshgrid(pos_w_x, pos_w_y, pos_w_z);
% % 
% %     % Interpolate coordinates on the contours of the wing and the body:
% %     
% %     x_sect_1_w = interp1(y_cont_w(1:((end+1)/2)),x_cont_w(1:((end+1)/2)),pos_w_y);
% %     x_sect_2_w = interp1(y_cont_w(((end+1)/2):end),x_cont_w(((end+1)/2):end),pos_w_y);
% %     
% %     
% %     [x_min, N_x_min] = min(x_cont_w);
% %     [x_max, N_x_max] = max(x_cont_w);
% %     
% %     N_x = length(x_cont_w);
% %     
% %     y_sect_1_w = interp1(x_cont_w(N_x_max:N_x_min),y_cont_w(N_x_max:N_x_min),pos_w_x);
% %     y_sect_2_w = interp1(x_cont_w([N_x:-1:N_x_min 1:N_x_max]),y_cont_w([N_x:-1:N_x_min 1:N_x_max]),pos_w_x);
% %     
% %     clear x_min x_max N_x_Min N_x_max
% %     
% %     Zw2 = zeros(sect_nr_w_y,sect_nr_w_x);
% %     
% %     
% %     for j = 1:sect_nr_w_y
% %         for i = 1:sect_nr_w_x
% %             if x_sect_1_w(j) >= Xw(j,i) && Xw(j,i) >= x_sect_2_w(j)
% %                 if y_sect_1_w(i) >= Yw(j,i) && Yw(j,i) >= y_sect_2_w(i)
% %                     Zw2(j,i) = 1;
% %                 end
% %             end
% %         end
% %     end
% %     
% %     Zh = zeros(sect_nr_w_y,sect_nr_w_x);
% %     
% %     
% %     for j = 1:sect_nr_w_y
% %         for i = 1:sect_nr_w_x
% %             if Zw2(j,i) == 1
% %                 
% %                 Zh(j,i) = sqrt(((x_sect_1_w(j)-x_sect_2_w(j))/2)^2-(pos_w_x(i)-(x_sect_1_w(j)-((x_sect_1_w(j)-x_sect_2_w(j))/2)))^2);
% %                 
% %             end
% %         end
% %     end
% %     
% %     
% % %     figure()
% % %     hold on
% % %     surf(Xw,Yw,Zh,'facecolor','none','edgecolor','k')
% % %     surf(Xw,Yw,-Zh,'facecolor','none','edgecolor','k')
% % %     axis equal
% % %     hold off
% % 
% % 
% % figure()
% % hold on
% % surf(J_pos_R(2)+wing_model.y_mod_L,-J_pos_R(1)+wing_model.x_mod_L,-J_pos_R(3)+wing_model.z_mod_L,'facecolor',[0.3 0.3 0.3],'edgecolor','none')
% % hold on
% % surf(-Yw,Xw,Zh,'facecolor','none','edgecolor','k')
% % surf(-Yw,Xw,-Zh,'facecolor','none','edgecolor','k')
% % for k = 1:sect_nr_w_x
% % plot3(-Yw2(:,:,k),Xw2(:,:,k),Zw(:,:,k),'o','MarkerSize',1.5,'MarkerEdgeColor', 'b')
% % end
% % hold off
% % title('Wing raster 3D')
% % xlabel('y [mm]')
% % ylabel('x [mm]')
% % zlabel('z [mm]')
% % xlim([-3.1 0])
% % ylim([-0.9 0.5])
% % zlim([-0.6 0.6])
% % axis equal
% % hold off

