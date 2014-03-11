function [ Fx_g, Fy_g, Fz_g, Mx_g, My_g, Mz_g ] = g_calib( theta_L, eta_L, phi_L, theta_R, eta_R, phi_R , g_grid_on )


    % Program to compute the forces and moments due to gravity in robofly.
    % First: load the calibration data and establish the calibration grid.
    % Second: use trilinear interpolation to establish a local value for
    % the gravitational forces and moments for a certain set of wing
    % kinematics.
    
    Fx_g = [];
    Fy_g = [];
    Fz_g = [];
    Mx_g = [];
    My_g = [];
    Mz_g = [];
    
    
    
    if g_grid_on == 0
    
        cd('C:/Users/Johan/Documents/Thesis_data/gravity_calibration/g_cal_07192013')
        
        temp_dir = dir;
    
        dir_names = {temp_dir.name};

        temp_isdir = [temp_dir.isdir];
        
        left_calib_names = [];
        
        right_calib_names = [];

        for j = 1:length(dir_names)

            if temp_isdir(j) == 0
                
                if isempty(findstr(char(dir_names(j)),'left')) == 0

                    left_calib_names = [left_calib_names; dir_names(j)];
                    
                elseif isempty(findstr(char(dir_names(j)),'right')) == 0
                    
                    right_calib_names = [right_calib_names; dir_names(j)];

                end

            end
        
        end
        
        nr_grid_pts = length(left_calib_names);
        
        left_calib_grid = zeros(9,nr_grid_pts);
        right_calib_grid = zeros(9,nr_grid_pts);
        
%         for j = 1:nr_grid_pts
%             
%             temp_left = load(char(left_calib_names(j)));
%             
%             left_calib_grid(1,j) = temp_left.deviation;
%             left_calib_grid(2,j) = temp_left.rotation;
%             left_calib_grid(3,j) = temp_left.stroke;
%             left_calib_grid(4,j) = 2.613e-5*mean(temp_left.ft(24000:end,3));
%             left_calib_grid(5,j) = 2.613e-5*mean(temp_left.ft(24000:end,1));
%             left_calib_grid(6,j) = 2.613e-5*mean(temp_left.ft(24000:end,5));
%             left_calib_grid(7,j) = 1e-3*2.613e-5*mean(temp_left.ft(24000:end,2));
%             left_calib_grid(8,j) = 1e-3*2.613e-5*mean(temp_left.ft(24000:end,4));
%             left_calib_grid(9,j) = 1e-3*2.613e-5*mean(temp_left.ft(24000:end,6));
%             
%             temp_right = load(char(right_calib_names(j)));
%             
%             right_calib_grid(1,j) = temp_right.deviation;
%             right_calib_grid(2,j) = temp_right.rotation;
%             right_calib_grid(3,j) = temp_right.stroke;
%             right_calib_grid(4,j) = 2.613e-5*mean(temp_right.ft(24000:end,3));
%             right_calib_grid(5,j) = 2.613e-5*mean(temp_right.ft(24000:end,1));
%             right_calib_grid(6,j) = 2.613e-5*mean(temp_right.ft(24000:end,5));
%             right_calib_grid(7,j) = 1e-3*2.613e-5*mean(temp_right.ft(24000:end,2));
%             right_calib_grid(8,j) = 1e-3*2.613e-5*mean(temp_right.ft(24000:end,4));
%             right_calib_grid(9,j) = 1e-3*2.613e-5*mean(temp_right.ft(24000:end,6));
%             
%         end
        
        for j = 1:nr_grid_pts
            
            temp_left = load(char(left_calib_names(j)));
            
            left_calib_grid(1,j) = temp_left.deviation;
            left_calib_grid(2,j) = temp_left.rotation;
            left_calib_grid(3,j) = temp_left.stroke;
            left_calib_grid(4,j) = mean(temp_left.ft(24000:end,3));
            left_calib_grid(5,j) = mean(temp_left.ft(24000:end,1));
            left_calib_grid(6,j) = mean(temp_left.ft(24000:end,5));
            left_calib_grid(7,j) = mean(temp_left.ft(24000:end,2));
            left_calib_grid(8,j) = mean(temp_left.ft(24000:end,4));
            left_calib_grid(9,j) = mean(temp_left.ft(24000:end,6));
            
            temp_right = load(char(right_calib_names(j)));
            
            right_calib_grid(1,j) = temp_right.deviation;
            right_calib_grid(2,j) = temp_right.rotation;
            right_calib_grid(3,j) = temp_right.stroke;
            right_calib_grid(4,j) = mean(temp_right.ft(24000:end,3));
            right_calib_grid(5,j) = mean(temp_right.ft(24000:end,1));
            right_calib_grid(6,j) = mean(temp_right.ft(24000:end,5));
            right_calib_grid(7,j) = mean(temp_right.ft(24000:end,2));
            right_calib_grid(8,j) = mean(temp_right.ft(24000:end,4));
            right_calib_grid(9,j) = mean(temp_right.ft(24000:end,6));
            
        end
        
%         figure()
%         plot3(left_calib_grid(3,:),left_calib_grid(1,:),left_calib_grid(2,:),'o')
        
        cd ..
        
        save('left_calib_grid.mat','left_calib_grid')
        
        save('right_calib_grid.mat','right_calib_grid')
                
    elseif g_grid_on == 1

                
        cd('C:/Users/Johan/Documents/Thesis_data/gravity_calibration')
        
               
        temp_L = load('left_calib_grid.mat');
        
        temp_R = load('right_calib_grid.mat');
        
        left_calib_grid = temp_L.left_calib_grid;
        
        right_calib_grid = temp_R.right_calib_grid;
        
        nr_grid_pts = length(left_calib_grid(1,:));
        
        theta_grid_left = unique(left_calib_grid(1,:));
        eta_grid_left = unique(left_calib_grid(2,:));
        phi_grid_left = unique(left_calib_grid(3,:));
        
        left_grid_Fx = zeros(length(theta_grid_left),length(eta_grid_left),length(phi_grid_left));
        left_grid_Fy = zeros(length(theta_grid_left),length(eta_grid_left),length(phi_grid_left));
        left_grid_Fz = zeros(length(theta_grid_left),length(eta_grid_left),length(phi_grid_left));
        left_grid_Mx = zeros(length(theta_grid_left),length(eta_grid_left),length(phi_grid_left));
        left_grid_My = zeros(length(theta_grid_left),length(eta_grid_left),length(phi_grid_left));
        left_grid_Mz = zeros(length(theta_grid_left),length(eta_grid_left),length(phi_grid_left));
        
        for j = 1:nr_grid_pts
            
            theta_temp = left_calib_grid(1,j);
            eta_temp = left_calib_grid(2,j);
            phi_temp = left_calib_grid(3,j);
            
            loc1 = find(theta_grid_left == theta_temp);
            loc2 = find(eta_grid_left == eta_temp);
            loc3 = find(phi_grid_left == phi_temp);
            
            left_grid_Fx(loc1,loc2,loc3) = left_calib_grid(4,j);
            left_grid_Fy(loc1,loc2,loc3) = left_calib_grid(5,j);
            left_grid_Fz(loc1,loc2,loc3) = left_calib_grid(6,j);
            left_grid_Mx(loc1,loc2,loc3) = left_calib_grid(7,j);
            left_grid_My(loc1,loc2,loc3) = left_calib_grid(8,j);
            left_grid_Mz(loc1,loc2,loc3) = left_calib_grid(9,j);
            
            
        end
        
        
        theta_grid_right = unique(right_calib_grid(1,:));
        eta_grid_right = unique(right_calib_grid(2,:));
        phi_grid_right = unique(right_calib_grid(3,:));
        
        right_grid_Fx = zeros(length(theta_grid_right),length(eta_grid_right),length(phi_grid_right));
        right_grid_Fy = zeros(length(theta_grid_right),length(eta_grid_right),length(phi_grid_right));
        right_grid_Fz = zeros(length(theta_grid_right),length(eta_grid_right),length(phi_grid_right));
        right_grid_Mx = zeros(length(theta_grid_right),length(eta_grid_right),length(phi_grid_right));
        right_grid_My = zeros(length(theta_grid_right),length(eta_grid_right),length(phi_grid_right));
        right_grid_Mz = zeros(length(theta_grid_right),length(eta_grid_right),length(phi_grid_right));
        
        for j = 1:nr_grid_pts
            
            theta_temp = right_calib_grid(1,j);
            eta_temp = right_calib_grid(2,j);
            phi_temp = right_calib_grid(3,j);
            
            loc1 = find(theta_grid_right == theta_temp);
            loc2 = find(eta_grid_right == eta_temp);
            loc3 = find(phi_grid_right == phi_temp);
            
            right_grid_Fx(loc1,loc2,loc3) = right_calib_grid(4,j);
            right_grid_Fy(loc1,loc2,loc3) = right_calib_grid(5,j);
            right_grid_Fz(loc1,loc2,loc3) = right_calib_grid(6,j);
            right_grid_Mx(loc1,loc2,loc3) = right_calib_grid(7,j);
            right_grid_My(loc1,loc2,loc3) = right_calib_grid(8,j);
            right_grid_Mz(loc1,loc2,loc3) = right_calib_grid(9,j);
            
        end          
        
        [eta_mesh_L,theta_mesh_L,phi_mesh_L] = meshgrid(eta_grid_left,theta_grid_left,phi_grid_left);
        
        [eta_mesh_R,theta_mesh_R,phi_mesh_R] = meshgrid(eta_grid_right,theta_grid_right,phi_grid_right);

        
%         Fx_g_left = interp3(theta_grid_left,eta_grid_left,phi_grid_left,left_grid_Fx,theta_grid_left,eta_grid_left,phi_grid_left);
%         Fy_g_left = interp3(theta_mesh_L,eta_mesh_L,phi_mesh_L,left_grid_Fy,theta_grid_left,eta_grid_left,phi_grid_left);
%         Fz_g_left = interp3(theta_mesh_L,eta_mesh_L,phi_mesh_L,left_grid_Fz,theta_L,eta_L,phi_L);
%         Mx_g_left = interp3(theta_mesh_L,eta_mesh_L,phi_mesh_L,left_grid_Mx,theta_L,eta_L,phi_L);
%         My_g_left = interp3(theta_mesh_L,eta_mesh_L,phi_mesh_L,left_grid_My,theta_L,eta_L,phi_L);
%         Mz_g_left = interp3(theta_mesh_L,eta_mesh_L,phi_mesh_L,left_grid_Mz,theta_L,eta_L,phi_L);
%         
%         Fx_g_right = interp3(theta_mesh_R,eta_mesh_R,phi_mesh_R,right_grid_Fx,theta_R,eta_L,phi_R);
%         Fy_g_right = interp3(theta_mesh_R,eta_mesh_R,phi_mesh_R,right_grid_Fy,theta_R,eta_L,phi_R);
%         Fz_g_right = interp3(theta_mesh_R,eta_mesh_R,phi_mesh_R,right_grid_Fz,theta_R,eta_L,phi_R);
%         Mx_g_right = interp3(theta_mesh_R,eta_mesh_R,phi_mesh_R,right_grid_Mx,theta_R,eta_L,phi_R);
%         My_g_right = interp3(theta_mesh_R,eta_mesh_R,phi_mesh_R,right_grid_My,theta_R,eta_L,phi_R);
%         Mz_g_right = interp3(theta_mesh_R,eta_mesh_R,phi_mesh_R,right_grid_Mz,theta_R,eta_L,phi_R);

        
        
        Fx_g_left = interp3(eta_mesh_L,theta_mesh_L,phi_mesh_L,left_grid_Fx,eta_L,theta_L,phi_L,'spline');
        Fy_g_left = interp3(eta_mesh_L,theta_mesh_L,phi_mesh_L,left_grid_Fy,eta_L,theta_L,phi_L,'spline');
        Fz_g_left = interp3(eta_mesh_L,theta_mesh_L,phi_mesh_L,left_grid_Fz,eta_L,theta_L,phi_L,'spline');
        Mx_g_left = interp3(eta_mesh_L,theta_mesh_L,phi_mesh_L,left_grid_Mx,eta_L,theta_L,phi_L,'spline');
        My_g_left = interp3(eta_mesh_L,theta_mesh_L,phi_mesh_L,left_grid_My,eta_L,theta_L,phi_L,'spline');
        Mz_g_left = interp3(eta_mesh_L,theta_mesh_L,phi_mesh_L,left_grid_Mz,eta_L,theta_L,phi_L,'spline');
        
        Fx_g_right = interp3(eta_mesh_R,theta_mesh_R,phi_mesh_R,right_grid_Fx,eta_R,theta_R,phi_R,'spline');
        Fy_g_right = interp3(eta_mesh_R,theta_mesh_R,phi_mesh_R,right_grid_Fy,eta_R,theta_R,phi_R,'spline');
        Fz_g_right = interp3(eta_mesh_R,theta_mesh_R,phi_mesh_R,right_grid_Fz,eta_R,theta_R,phi_R,'spline');
        Mx_g_right = interp3(eta_mesh_R,theta_mesh_R,phi_mesh_R,right_grid_Mx,eta_R,theta_R,phi_R,'spline');
        My_g_right = interp3(eta_mesh_R,theta_mesh_R,phi_mesh_R,right_grid_My,eta_R,theta_R,phi_R,'spline');
        Mz_g_right = interp3(eta_mesh_R,theta_mesh_R,phi_mesh_R,right_grid_Mz,eta_R,theta_R,phi_R,'spline');
        
        
        Fx_g = Fx_g_left + Fx_g_right;
        Fy_g = Fy_g_left + Fy_g_right;
        Fz_g = Fz_g_left + Fz_g_right;
        Mx_g = Mx_g_left + Mx_g_right;
        My_g = My_g_left + My_g_right;
        Mz_g = Mz_g_left + Mz_g_right;
        

        
    end
        
        
    
    

end

