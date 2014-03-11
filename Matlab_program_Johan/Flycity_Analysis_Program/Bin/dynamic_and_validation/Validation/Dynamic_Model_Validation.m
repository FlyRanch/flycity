function Dynamic_Model_Validation( settings, pathDB )

    savefile1 = 'pathDB8.mat';
    
%     savefile2 = 'DynSim_temp.mat';

    savefile2 = 'DynSim_temp.mat';

%     savefile2 = 'DynSim_temp_FMI_off.mat';

    nr_seq = settings.nr_of_seq;
    
    K = 25; % Number of steady flight wingbeats
    
    N = 75; % Number of random maneuvering flight wingbeats
    
    % Create a list of names containing maneuvering wingbeats:
    
    ax_names = fieldnames(pathDB.maneuver.ax);
    ay_names = fieldnames(pathDB.maneuver.ay);
    az_names = fieldnames(pathDB.maneuver.az);
    wx_names = fieldnames(pathDB.maneuver.wx);
    wy_names = fieldnames(pathDB.maneuver.wy);
    wz_names = fieldnames(pathDB.maneuver.wz);
    
    
    
    
    % List maneuvering wingbeats:
    
    man_wb_names = {};  
    
    for k = 1:length(ax_names)
        
        man_wb_names{end+1} = ['ax_' char(ax_names(k))];
        
    end
    
    for k = 1:length(ay_names)
        
        man_wb_names{end+1} = ['ay_' char(ay_names(k))];
        
    end
    
    for k = 1:length(az_names)
        
        man_wb_names{end+1} = ['az_' char(az_names(k))];
        
    end
    
    for k = 1:length(wx_names)
        
        man_wb_names{end+1} = ['wx_' char(wx_names(k))];
        
    end 
    
    for k = 1:length(wy_names)
        
        man_wb_names{end+1} = ['wy_' char(wy_names(k))];
        
    end 
    
    for k = 1:length(wz_names)
        
        man_wb_names{end+1} = ['wz_' char(wz_names(k))];
        
    end 
    
    % Test if simulation data is available:
    
    if exist(char(savefile2))==2
        
        temp = load(char(savefile2));
        
        if isfield(temp,'DynSim_avg')==1
        
            DynSim_avg = temp.DynSim_avg;
        
        end
        
        if isfield(temp,'DynSim_man')==1
            
            DynSim_man = temp.DynSim_man;
            
        end
        
        if isfield(temp,'DynSim_steady')==1
            
            DynSim_steady = temp.DynSim_steady;
            
        end
        
        clear temp
        
    end
    
    % Create a list of sequence numbers and wb numbers of maneuvering
    % wingbeats:
    
    r = length(man_wb_names);
    
    man_seq_wb = zeros(r,2);
    
    for j = 1:r
        
         if strncmpi(man_wb_names(j),'ax',2) == 1
             
            t_string    = char(man_wb_names(j));
            wb_nr_t2    = str2num(t_string(7:end)); 
            wb_nr_t     = pathDB.maneuver.ax.(char(['wb_' num2str(wb_nr_t2)])).wb_nr;
            seq_nr_t    = pathDB.maneuver.ax.(char(['wb_' num2str(wb_nr_t2)])).seq_nr;
            
        end
                    
        if strncmpi(man_wb_names(j),'ay',2) == 1
            
            t_string    = char(man_wb_names(j));
            wb_nr_t2    = str2num(t_string(7:end)); 
            wb_nr_t     = pathDB.maneuver.ay.(char(['wb_' num2str(wb_nr_t2)])).wb_nr;
            seq_nr_t    = pathDB.maneuver.ay.(char(['wb_' num2str(wb_nr_t2)])).seq_nr;
                        
        end
        
        if strncmpi(man_wb_names(j),'az',2) == 1
            
            t_string    = char(man_wb_names(j));
            wb_nr_t2    = str2num(t_string(7:end)); 
            wb_nr_t     = pathDB.maneuver.az.(char(['wb_' num2str(wb_nr_t2)])).wb_nr;
            seq_nr_t    = pathDB.maneuver.az.(char(['wb_' num2str(wb_nr_t2)])).seq_nr;
            
        end
                    
        if strncmpi(man_wb_names(j),'wx',2) == 1
            
            t_string    = char(man_wb_names(j));
            wb_nr_t2    = str2num(t_string(7:end)); 
            wb_nr_t     = pathDB.maneuver.wx.(char(['wb_' num2str(wb_nr_t2)])).wb_nr;
            seq_nr_t    = pathDB.maneuver.wx.(char(['wb_' num2str(wb_nr_t2)])).seq_nr;
                        
        end
        
        if strncmpi(man_wb_names(j),'wy',2) == 1
            
            t_string    = char(man_wb_names(j));
            wb_nr_t2    = str2num(t_string(7:end)); 
            wb_nr_t     = pathDB.maneuver.wy.(char(['wb_' num2str(wb_nr_t2)])).wb_nr;
            seq_nr_t    = pathDB.maneuver.wy.(char(['wb_' num2str(wb_nr_t2)])).seq_nr;
                        
        end
        
        if strncmpi(man_wb_names(j),'wz',2) == 1
            
            t_string    = char(man_wb_names(j));
            wb_nr_t2    = str2num(t_string(7:end)); 
            wb_nr_t     = pathDB.maneuver.wz.(char(['wb_' num2str(wb_nr_t2)])).wb_nr;
            seq_nr_t    = pathDB.maneuver.wz.(char(['wb_' num2str(wb_nr_t2)])).seq_nr;
                        
        end
        
        man_seq_wb(j,1) = wb_nr_t;
        man_seq_wb(j,2) = seq_nr_t;
        
    end
    
        
    if exist('DynSim_avg')==1
        
        'DynSim_avg exists'
                
    else

        'Creating DynSim_avg'

        % Perform the Dynamic Simulation on all the average wingbeats:

        DynSim_avg = {};

        % All average wingbeats

        for i = 1:nr_seq

            if isnan(pathDB.poly_fit.a_avg.f(i))==0

                i

                [kine_avg, sim_data_avg]    = Dynamic_simulation_avg( settings, pathDB, i, 0, 1);

                DynSim_avg.(char(['seq_' num2str(i)])).kine         = kine_avg;
                DynSim_avg.(char(['seq_' num2str(i)])).sim_data     = sim_data_avg;

            end

        end
        
        save(savefile2,'DynSim_avg')
        
    
    end
    

    
    % Obtain K steady wingbeats:
    
    steady_wb_names = fieldnames(pathDB.rand_wbs);
    
    q = length(steady_wb_names);
    
    if exist('DynSim_steady')==1
        
        'DynSim_steady exists'
        
        DynSim_steady_names = fieldnames(DynSim_steady);
        
        Q = length(DynSim_steady_names);
        
        R = K-Q;
        
        i = R;
        
        if R > 0
            
            while i >= 0
                
                i
            
                % Pick a random number:
                
                j = round(rand*q);
                
                if j == 0
                    
                    j = round(rand*q);
                    
                end
                
                j
                
                if sum(strcmp(DynSim_steady_names,char(steady_wb_names(j)))) == 0
                    
                    wb_nr_steady  = pathDB.rand_wbs.(char(['wb_' num2str(j)])).wb_nr;
                    seq_nr_steady = pathDB.rand_wbs.(char(['wb_' num2str(j)])).seq_nr;
                    
                    is_man = 0;
                    
                    if isfield(DynSim_avg , (char(['seq_' char(seq_nr_steady)]))) == 1
                    
                        for k = 1:size(man_seq_wb,1)

                            if man_seq_wb(k,1)==wb_nr_steady && man_seq_wb(k,2)==seq_nr_steady

                                is_man = 1;

                            end

                        end

                        if is_man == 0


                            body_IC.xyz_0   = [0; 0; 0];
                            body_IC.qb_0    = pathDB.rand_wbs.(char(['wb_' num2str(j)])).body_kin.qb(1,:)';
                            body_IC.vb_0    = quat2mat(body_IC.qb_0)*pathDB.rand_wbs.(char(['wb_' num2str(j)])).body_kin.uvw(1,:)';
                            body_IC.wb_0    = pathDB.rand_wbs.(char(['wb_' num2str(j)])).body_kin.wb(1,:)';
                            body_IC.F_0     = DynSim_avg.(char(['seq_' num2str(seq_nr_steady)])).sim_data.F_0;
                            body_IC.M_0     = DynSim_avg.(char(['seq_' num2str(seq_nr_steady)])).sim_data.M_0;
%                             body_IC.F_0     = [0; 0; 0];
%                             body_IC.M_0     = [0; 0; 0];

                            [kine_steady, sim_data_steady] = Dynamic_simulation_steady( settings, pathDB, body_IC, char(['wb_' num2str(j)]), seq_nr_steady, 1 );

                            DynSim_steady.(char(['wb_' num2str(j)])).kine         = kine_steady;
                            DynSim_steady.(char(['wb_' num2str(j)])).sim_data     = sim_data_steady;

                            i = i-1;

                        end
                    
                    end
                    
                end
                
            end
            
        end
        
        save(savefile2,'DynSim_avg','DynSim_steady')
        
    else
   
        'Creating DynSim_steady'
        
        i = 0;
        
        while i < K
            
            i
            
            % Pick a random number:
            
            j = round(rand*q);
            
            if j == 0
                    
                j = round(rand*q);
                    
            end
                
            j
            
            wb_nr_steady  = pathDB.rand_wbs.(char(['wb_' num2str(j)])).wb_nr;
            seq_nr_steady = pathDB.rand_wbs.(char(['wb_' num2str(j)])).seq_nr;

            is_man = 0;
            
            if isfield(DynSim_avg , (char(['seq_' char(seq_nr_steady)]))) == 1
                    
                for k = 1:size(man_seq_wb,1)

                    if man_seq_wb(k,1)==wb_nr_steady && man_seq_wb(k,2)==seq_nr_steady

                        is_man = 1;

                    end

                end

                if is_man == 0
                                                       
                    body_IC.xyz_0   = [0; 0; 0];
                    body_IC.qb_0    = pathDB.rand_wbs.(char(['wb_' num2str(j)])).body_kin.qb(1,:)';
                    body_IC.vb_0    = quat2mat(body_IC.qb_0)*pathDB.rand_wbs.(char(['wb_' num2str(j)])).body_kin.uvw(1,:)';
                    body_IC.wb_0    = pathDB.rand_wbs.(char(['wb_' num2str(j)])).body_kin.wb(1,:)';
                    body_IC.F_0     = DynSim_avg.(char(['seq_' num2str(seq_nr_steady)])).sim_data.F_0;
                    body_IC.M_0     = DynSim_avg.(char(['seq_' num2str(seq_nr_steady)])).sim_data.M_0;
%                     body_IC.F_0     = [0; 0; 0];
%                     body_IC.M_0     = [0; 0; 0];                    

                    [kine_steady, sim_data_steady] = Dynamic_simulation_steady( settings, pathDB, body_IC, char(['wb_' num2str(j)]), seq_nr_steady, 1 );

                    DynSim_steady.(char(['wb_' num2str(j)])).kine         = kine_steady;
                    DynSim_steady.(char(['wb_' num2str(j)])).sim_data     = sim_data_steady;

                    i = i+1;
                    
                end  
            
            end
            
        end
        
        save(savefile2,'DynSim_avg','DynSim_steady')
        
    end
    
    % Perform the Dynamic Simulation on N maneuvering wingbeats:
    
            
    if exist('DynSim_man')==1
        
        'DynSim_man exists'
        
        DynSim_man_names = fieldnames(DynSim_man);
        
        P = length(DynSim_man_names);
        
        M = N-P;
        
        i = M;
        
        if M > 0
            
            while i >= 0
                
                i
                
                % Pick a random number:
                
                j = round(rand*r);
                
                if j == 0
                    
                    j = round(rand*r);
                    
                end
                
                j
                
                if sum(strcmp(DynSim_man_names,char(man_wb_names(j)))) == 0
                    
                    if strncmpi(man_wb_names(j),'ax',2) == 1
                        
                        t_string = char(man_wb_names(j));
                        wb_nr = str2num(t_string(7:end));
                        
                        seq_nr          = pathDB.maneuver.ax.(char(['wb_' num2str(wb_nr)])).seq_nr;
                        body_IC.xyz_0   = [0; 0; 0];
                        body_IC.qb_0    = pathDB.maneuver.ax.(char(['wb_' num2str(wb_nr)])).body_kin.qb(1,:)';
                        body_IC.vb_0    = quat2mat(body_IC.qb_0)*pathDB.maneuver.ax.(char(['wb_' num2str(wb_nr)])).body_kin.uvw(1,:)';
                        body_IC.wb_0    = pathDB.maneuver.ax.(char(['wb_' num2str(wb_nr)])).body_kin.wb(1,:)';
                        body_IC.F_0     = DynSim_avg.(char(['seq_' num2str(seq_nr)])).sim_data.F_0;
                        body_IC.M_0     = DynSim_avg.(char(['seq_' num2str(seq_nr)])).sim_data.M_0;
%                         body_IC.F_0     = [0; 0; 0];
%                         body_IC.M_0     = [0; 0; 0];                        
                        
                        [kine_man, sim_data_man] = Dynamic_simulation_man( settings, pathDB, body_IC, 'ax', char(ax_names(wb_nr)), seq_nr, 1 );
                        
                        DynSim_man.(char(man_wb_names(j))).kine         = kine_man;
                        DynSim_man.(char(man_wb_names(j))).sim_data     = sim_data_man;
                        
                    end
                    
                    if strncmpi(man_wb_names(j),'ay',2) == 1
                        
                        t_string = char(man_wb_names(j));
                        wb_nr = str2num(t_string(7:end));
                        
                        seq_nr          = pathDB.maneuver.ay.(char(['wb_' num2str(wb_nr)])).seq_nr;
                        body_IC.xyz_0   = [0; 0; 0];
                        body_IC.qb_0    = pathDB.maneuver.ay.(char(['wb_' num2str(wb_nr)])).body_kin.qb(1,:)';
                        body_IC.vb_0    = quat2mat(body_IC.qb_0)*pathDB.maneuver.ay.(char(['wb_' num2str(wb_nr)])).body_kin.uvw(1,:)';
                        body_IC.wb_0    = pathDB.maneuver.ay.(char(['wb_' num2str(wb_nr)])).body_kin.wb(1,:)';
                        body_IC.F_0     = DynSim_avg.(char(['seq_' num2str(seq_nr)])).sim_data.F_0;
                        body_IC.M_0     = DynSim_avg.(char(['seq_' num2str(seq_nr)])).sim_data.M_0;
%                         body_IC.F_0     = [0; 0; 0];
%                         body_IC.M_0     = [0; 0; 0];                          
                        
                        [kine_man, sim_data_man] = Dynamic_simulation_man( settings, pathDB, body_IC, 'ay', char(ay_names(wb_nr)), seq_nr, 1 );
                        
                        DynSim_man.(char(man_wb_names(j))).kine         = kine_man;
                        DynSim_man.(char(man_wb_names(j))).sim_data     = sim_data_man;
                        
                    end
                    
                    if strncmpi(man_wb_names(j),'az',2) == 1
                        
                        t_string = char(man_wb_names(j));
                        wb_nr = str2num(t_string(7:end));
                        
                        seq_nr          = pathDB.maneuver.az.(char(['wb_' num2str(wb_nr)])).seq_nr;
                        body_IC.xyz_0   = [0; 0; 0];
                        body_IC.qb_0    = pathDB.maneuver.az.(char(['wb_' num2str(wb_nr)])).body_kin.qb(1,:)';
                        body_IC.vb_0    = quat2mat(body_IC.qb_0)*pathDB.maneuver.az.(char(['wb_' num2str(wb_nr)])).body_kin.uvw(1,:)';
                        body_IC.wb_0    = pathDB.maneuver.az.(char(['wb_' num2str(wb_nr)])).body_kin.wb(1,:)';
                        body_IC.F_0     = DynSim_avg.(char(['seq_' num2str(seq_nr)])).sim_data.F_0;
                        body_IC.M_0     = DynSim_avg.(char(['seq_' num2str(seq_nr)])).sim_data.M_0;
%                         body_IC.F_0     = [0; 0; 0];
%                         body_IC.M_0     = [0; 0; 0];                          
                        
                        [kine_man, sim_data_man] = Dynamic_simulation_man( settings, pathDB, body_IC, 'az', char(az_names(wb_nr)), seq_nr, 1 );
                        
                        DynSim_man.(char(man_wb_names(j))).kine         = kine_man;
                        DynSim_man.(char(man_wb_names(j))).sim_data     = sim_data_man;
                        
                    end
                    
                    if strncmpi(man_wb_names(j),'wx',2) == 1
                        
                        t_string = char(man_wb_names(j));
                        wb_nr = str2num(t_string(7:end));
                        
                        seq_nr          = pathDB.maneuver.wx.(char(['wb_' num2str(wb_nr)])).seq_nr;
                        body_IC.xyz_0   = [0; 0; 0];
                        body_IC.qb_0    = pathDB.maneuver.wx.(char(['wb_' num2str(wb_nr)])).body_kin.qb(1,:)';
                        body_IC.vb_0    = quat2mat(body_IC.qb_0)*pathDB.maneuver.wx.(char(['wb_' num2str(wb_nr)])).body_kin.uvw(1,:)';
                        body_IC.wb_0    = pathDB.maneuver.wx.(char(['wb_' num2str(wb_nr)])).body_kin.wb(1,:)';
                        body_IC.F_0     = DynSim_avg.(char(['seq_' num2str(seq_nr)])).sim_data.F_0;
                        body_IC.M_0     = DynSim_avg.(char(['seq_' num2str(seq_nr)])).sim_data.M_0;
%                         body_IC.F_0     = [0; 0; 0];
%                         body_IC.M_0     = [0; 0; 0];                          
                        
                        [kine_man, sim_data_man] = Dynamic_simulation_man( settings, pathDB, body_IC, 'wx', char(wx_names(wb_nr)), seq_nr, 1 );
                        
                        DynSim_man.(char(man_wb_names(j))).kine         = kine_man;
                        DynSim_man.(char(man_wb_names(j))).sim_data     = sim_data_man;
                        
                    end
                    
                    if strncmpi(man_wb_names(j),'wy',2) == 1
                        
                        t_string = char(man_wb_names(j));
                        wb_nr = str2num(t_string(7:end));
                        
                        seq_nr          = pathDB.maneuver.wy.(char(['wb_' num2str(wb_nr)])).seq_nr;
                        body_IC.xyz_0   = [0; 0; 0];
                        body_IC.qb_0    = pathDB.maneuver.wy.(char(['wb_' num2str(wb_nr)])).body_kin.qb(1,:)';
                        body_IC.vb_0    = quat2mat(body_IC.qb_0)*pathDB.maneuver.wy.(char(['wb_' num2str(wb_nr)])).body_kin.uvw(1,:)';
                        body_IC.wb_0    = pathDB.maneuver.wy.(char(['wb_' num2str(wb_nr)])).body_kin.wb(1,:)';
                        body_IC.F_0     = DynSim_avg.(char(['seq_' num2str(seq_nr)])).sim_data.F_0;
                        body_IC.M_0     = DynSim_avg.(char(['seq_' num2str(seq_nr)])).sim_data.M_0;
%                         body_IC.F_0     = [0; 0; 0];
%                         body_IC.M_0     = [0; 0; 0];                          
                        
                        [kine_man, sim_data_man] = Dynamic_simulation_man( settings, pathDB, body_IC, 'wy', char(wy_names(wb_nr)), seq_nr, 1 );
                        
                        DynSim_man.(char(man_wb_names(j))).kine         = kine_man;
                        DynSim_man.(char(man_wb_names(j))).sim_data     = sim_data_man;
                        
                    end
                    
                    if strncmpi(man_wb_names(j),'wz',2) == 1
                        
                        t_string = char(man_wb_names(j));
                        wb_nr = str2num(t_string(7:end));
                        
                        seq_nr          = pathDB.maneuver.wz.(char(['wb_' num2str(wb_nr)])).seq_nr;
                        body_IC.xyz_0   = [0; 0; 0];
                        body_IC.qb_0    = pathDB.maneuver.wz.(char(['wb_' num2str(wb_nr)])).body_kin.qb(1,:)';
                        body_IC.vb_0    = quat2mat(body_IC.qb_0)*pathDB.maneuver.wz.(char(['wb_' num2str(wb_nr)])).body_kin.uvw(1,:)';
                        body_IC.wb_0    = pathDB.maneuver.wz.(char(['wb_' num2str(wb_nr)])).body_kin.wb(1,:)';
                        body_IC.F_0     = DynSim_avg.(char(['seq_' num2str(seq_nr)])).sim_data.F_0;
                        body_IC.M_0     = DynSim_avg.(char(['seq_' num2str(seq_nr)])).sim_data.M_0;
%                         body_IC.F_0     = [0; 0; 0];
%                         body_IC.M_0     = [0; 0; 0];                          
                        
                        [kine_man, sim_data_man] = Dynamic_simulation_man( settings, pathDB, body_IC, 'wz', char(wz_names(wb_nr)), seq_nr, 1 );
                        
                        DynSim_man.(char(man_wb_names(j))).kine         = kine_man;
                        DynSim_man.(char(man_wb_names(j))).sim_data     = sim_data_man;
                        
                    end
                    
                    i = i-1;
                    
                end
                
            end
            
        end
        
        save(savefile2,'DynSim_avg','DynSim_steady','DynSim_man')
        
    else
        
        'Creating DynSim_man'
        
        i = 0;
        
        while i < N
                
                i
                
                % Pick a random number:
                
                j = round(rand*r);
                
                if j == 0
                    
                    j = round(rand*r);
                    
                end
                
                j
                
                if i == 0
                                    
                    if strncmpi(man_wb_names(j),'ax',2) == 1
                        
                        t_string = char(man_wb_names(j));
                        wb_nr = str2num(t_string(7:end));
                        
                        seq_nr          = pathDB.maneuver.ax.(char(['wb_' num2str(wb_nr)])).seq_nr;
                        body_IC.xyz_0   = [0; 0; 0];
                        body_IC.qb_0    = pathDB.maneuver.ax.(char(['wb_' num2str(wb_nr)])).body_kin.qb(1,:)';
                        body_IC.vb_0    = quat2mat(body_IC.qb_0)*pathDB.maneuver.ax.(char(['wb_' num2str(wb_nr)])).body_kin.uvw(1,:)';
                        body_IC.wb_0    = pathDB.maneuver.ax.(char(['wb_' num2str(wb_nr)])).body_kin.wb(1,:)';
                        body_IC.F_0     = DynSim_avg.(char(['seq_' num2str(seq_nr)])).sim_data.F_0;
                        body_IC.M_0     = DynSim_avg.(char(['seq_' num2str(seq_nr)])).sim_data.M_0;
%                         body_IC.F_0     = [0; 0; 0];
%                         body_IC.M_0     = [0; 0; 0];                          
                        
                        [kine_man, sim_data_man] = Dynamic_simulation_man( settings, pathDB, body_IC, 'ax', char(ax_names(wb_nr)), seq_nr, 1 );
                        
                        DynSim_man.(char(man_wb_names(j))).kine         = kine_man;
                        DynSim_man.(char(man_wb_names(j))).sim_data     = sim_data_man;
                        
                    end
                    
                    if strncmpi(man_wb_names(j),'ay',2) == 1
                        
                        t_string = char(man_wb_names(j));
                        wb_nr = str2num(t_string(7:end));
                        
                        seq_nr          = pathDB.maneuver.ay.(char(['wb_' num2str(wb_nr)])).seq_nr;
                        body_IC.xyz_0   = [0; 0; 0];
                        body_IC.qb_0    = pathDB.maneuver.ay.(char(['wb_' num2str(wb_nr)])).body_kin.qb(1,:)';
                        body_IC.vb_0    = quat2mat(body_IC.qb_0)*pathDB.maneuver.ay.(char(['wb_' num2str(wb_nr)])).body_kin.uvw(1,:)';
                        body_IC.wb_0    = pathDB.maneuver.ay.(char(['wb_' num2str(wb_nr)])).body_kin.wb(1,:)';
                        body_IC.F_0     = DynSim_avg.(char(['seq_' num2str(seq_nr)])).sim_data.F_0;
                        body_IC.M_0     = DynSim_avg.(char(['seq_' num2str(seq_nr)])).sim_data.M_0;
%                         body_IC.F_0     = [0; 0; 0];
%                         body_IC.M_0     = [0; 0; 0];                          
                        
                        [kine_man, sim_data_man] = Dynamic_simulation_man( settings, pathDB, body_IC, 'ay', char(ay_names(wb_nr)), seq_nr, 1 );
                        
                        DynSim_man.(char(man_wb_names(j))).kine         = kine_man;
                        DynSim_man.(char(man_wb_names(j))).sim_data     = sim_data_man;
                        
                    end
                    
                    if strncmpi(man_wb_names(j),'az',2) == 1
                        
                        t_string = char(man_wb_names(j));
                        wb_nr = str2num(t_string(7:end));
                        
                        seq_nr          = pathDB.maneuver.az.(char(['wb_' num2str(wb_nr)])).seq_nr;
                        body_IC.xyz_0   = [0; 0; 0];
                        body_IC.qb_0    = pathDB.maneuver.az.(char(['wb_' num2str(wb_nr)])).body_kin.qb(1,:)';
                        body_IC.vb_0    = quat2mat(body_IC.qb_0)*pathDB.maneuver.az.(char(['wb_' num2str(wb_nr)])).body_kin.uvw(1,:)';
                        body_IC.wb_0    = pathDB.maneuver.az.(char(['wb_' num2str(wb_nr)])).body_kin.wb(1,:)';
                        body_IC.F_0     = DynSim_avg.(char(['seq_' num2str(seq_nr)])).sim_data.F_0;
                        body_IC.M_0     = DynSim_avg.(char(['seq_' num2str(seq_nr)])).sim_data.M_0;
%                         body_IC.F_0     = [0; 0; 0];
%                         body_IC.M_0     = [0; 0; 0];                          
                        
                        [kine_man, sim_data_man] = Dynamic_simulation_man( settings, pathDB, body_IC, 'az', char(az_names(wb_nr)), seq_nr, 1 );
                        
                        DynSim_man.(char(man_wb_names(j))).kine         = kine_man;
                        DynSim_man.(char(man_wb_names(j))).sim_data     = sim_data_man;
                        
                    end
                    
                    if strncmpi(man_wb_names(j),'wx',2) == 1
                        
                        t_string = char(man_wb_names(j));
                        wb_nr = str2num(t_string(7:end));
                        
                        seq_nr          = pathDB.maneuver.wx.(char(['wb_' num2str(wb_nr)])).seq_nr;
                        body_IC.xyz_0   = [0; 0; 0];
                        body_IC.qb_0    = pathDB.maneuver.wx.(char(['wb_' num2str(wb_nr)])).body_kin.qb(1,:)';
                        body_IC.vb_0    = quat2mat(body_IC.qb_0)*pathDB.maneuver.wx.(char(['wb_' num2str(wb_nr)])).body_kin.uvw(1,:)';
                        body_IC.wb_0    = pathDB.maneuver.wx.(char(['wb_' num2str(wb_nr)])).body_kin.wb(1,:)';
                        body_IC.F_0     = DynSim_avg.(char(['seq_' num2str(seq_nr)])).sim_data.F_0;
                        body_IC.M_0     = DynSim_avg.(char(['seq_' num2str(seq_nr)])).sim_data.M_0;
%                         body_IC.F_0     = [0; 0; 0];
%                         body_IC.M_0     = [0; 0; 0];                          
                        
                        [kine_man, sim_data_man] = Dynamic_simulation_man( settings, pathDB, body_IC, 'wx', char(wx_names(wb_nr)), seq_nr, 1 );
                        
                        DynSim_man.(char(man_wb_names(j))).kine         = kine_man;
                        DynSim_man.(char(man_wb_names(j))).sim_data     = sim_data_man;
                        
                    end
                    
                    if strncmpi(man_wb_names(j),'wy',2) == 1
                        
                        t_string = char(man_wb_names(j));
                        wb_nr = str2num(t_string(7:end));
                        
                        seq_nr          = pathDB.maneuver.wy.(char(['wb_' num2str(wb_nr)])).seq_nr;
                        body_IC.xyz_0   = [0; 0; 0];
                        body_IC.qb_0    = pathDB.maneuver.wy.(char(['wb_' num2str(wb_nr)])).body_kin.qb(1,:)';
                        body_IC.vb_0    = quat2mat(body_IC.qb_0)*pathDB.maneuver.wy.(char(['wb_' num2str(wb_nr)])).body_kin.uvw(1,:)';
                        body_IC.wb_0    = pathDB.maneuver.wy.(char(['wb_' num2str(wb_nr)])).body_kin.wb(1,:)';
                        body_IC.F_0     = DynSim_avg.(char(['seq_' num2str(seq_nr)])).sim_data.F_0;
                        body_IC.M_0     = DynSim_avg.(char(['seq_' num2str(seq_nr)])).sim_data.M_0;
%                         body_IC.F_0     = [0; 0; 0];
%                         body_IC.M_0     = [0; 0; 0];                          
                        
                        [kine_man, sim_data_man] = Dynamic_simulation_man( settings, pathDB, body_IC, 'wy', char(wy_names(wb_nr)), seq_nr, 1 );
                        
                        DynSim_man.(char(man_wb_names(j))).kine         = kine_man;
                        DynSim_man.(char(man_wb_names(j))).sim_data     = sim_data_man;
                        
                    end
                    
                    if strncmpi(man_wb_names(j),'wz',2) == 1
                        
                        t_string = char(man_wb_names(j));
                        wb_nr = str2num(t_string(7:end));
                        
                        seq_nr          = pathDB.maneuver.wz.(char(['wb_' num2str(wb_nr)])).seq_nr;
                        body_IC.xyz_0   = [0; 0; 0];
                        body_IC.qb_0    = pathDB.maneuver.wz.(char(['wb_' num2str(wb_nr)])).body_kin.qb(1,:)';
                        body_IC.vb_0    = quat2mat(body_IC.qb_0)*pathDB.maneuver.wz.(char(['wb_' num2str(wb_nr)])).body_kin.uvw(1,:)';
                        body_IC.wb_0    = pathDB.maneuver.wz.(char(['wb_' num2str(wb_nr)])).body_kin.wb(1,:)';
                        body_IC.F_0     = DynSim_avg.(char(['seq_' num2str(seq_nr)])).sim_data.F_0;
                        body_IC.M_0     = DynSim_avg.(char(['seq_' num2str(seq_nr)])).sim_data.M_0;
%                         body_IC.F_0     = [0; 0; 0];
%                         body_IC.M_0     = [0; 0; 0];                          
                        
                        [kine_man, sim_data_man] = Dynamic_simulation_man( settings, pathDB, body_IC, 'wz', char(wz_names(wb_nr)), seq_nr, 1 );
                        
                        DynSim_man.(char(man_wb_names(j))).kine         = kine_man;
                        DynSim_man.(char(man_wb_names(j))).sim_data     = sim_data_man;
                        
                    end
                    
                    i = i+1;
                    
                elseif sum(strcmp(fieldnames(DynSim_man),char(man_wb_names(j)))) == 0
                    
                    if strncmpi(man_wb_names(j),'ax',2) == 1
                        
                        t_string = char(man_wb_names(j));
                        wb_nr = str2num(t_string(7:end));
                        
                        seq_nr          = pathDB.maneuver.ax.(char(['wb_' num2str(wb_nr)])).seq_nr;
                        body_IC.xyz_0   = [0; 0; 0];
                        body_IC.qb_0    = pathDB.maneuver.ax.(char(['wb_' num2str(wb_nr)])).body_kin.qb(1,:)';
                        body_IC.vb_0    = quat2mat(body_IC.qb_0)*pathDB.maneuver.ax.(char(['wb_' num2str(wb_nr)])).body_kin.uvw(1,:)';
                        body_IC.wb_0    = pathDB.maneuver.ax.(char(['wb_' num2str(wb_nr)])).body_kin.wb(1,:)';
                        body_IC.F_0     = DynSim_avg.(char(['seq_' num2str(seq_nr)])).sim_data.F_0;
                        body_IC.M_0     = DynSim_avg.(char(['seq_' num2str(seq_nr)])).sim_data.M_0;
%                         body_IC.F_0     = [0; 0; 0];
%                         body_IC.M_0     = [0; 0; 0];                          
                        
                        [kine_man, sim_data_man] = Dynamic_simulation_man( settings, pathDB, body_IC, 'ax', char(ax_names(wb_nr)), seq_nr, 1 );
                        
                        DynSim_man.(char(man_wb_names(j))).kine         = kine_man;
                        DynSim_man.(char(man_wb_names(j))).sim_data     = sim_data_man;
                        
                    end
                    
                    if strncmpi(man_wb_names(j),'ay',2) == 1
                        
                        t_string = char(man_wb_names(j));
                        wb_nr = str2num(t_string(7:end));
                        
                        seq_nr          = pathDB.maneuver.ay.(char(['wb_' num2str(wb_nr)])).seq_nr;
                        body_IC.xyz_0   = [0; 0; 0];
                        body_IC.qb_0    = pathDB.maneuver.ay.(char(['wb_' num2str(wb_nr)])).body_kin.qb(1,:)';
                        body_IC.vb_0    = quat2mat(body_IC.qb_0)*pathDB.maneuver.ay.(char(['wb_' num2str(wb_nr)])).body_kin.uvw(1,:)';
                        body_IC.wb_0    = pathDB.maneuver.ay.(char(['wb_' num2str(wb_nr)])).body_kin.wb(1,:)';
                        body_IC.F_0     = DynSim_avg.(char(['seq_' num2str(seq_nr)])).sim_data.F_0;
                        body_IC.M_0     = DynSim_avg.(char(['seq_' num2str(seq_nr)])).sim_data.M_0;
%                         body_IC.F_0     = [0; 0; 0];
%                         body_IC.M_0     = [0; 0; 0];                          
                        
                        [kine_man, sim_data_man] = Dynamic_simulation_man( settings, pathDB, body_IC, 'ay', char(ay_names(wb_nr)), seq_nr, 1 );
                        
                        DynSim_man.(char(man_wb_names(j))).kine         = kine_man;
                        DynSim_man.(char(man_wb_names(j))).sim_data     = sim_data_man;
                        
                    end
                    
                    if strncmpi(man_wb_names(j),'az',2) == 1
                        
                        t_string = char(man_wb_names(j));
                        wb_nr = str2num(t_string(7:end));
                        
                        seq_nr          = pathDB.maneuver.az.(char(['wb_' num2str(wb_nr)])).seq_nr;
                        body_IC.xyz_0   = [0; 0; 0];
                        body_IC.qb_0    = pathDB.maneuver.az.(char(['wb_' num2str(wb_nr)])).body_kin.qb(1,:)';
                        body_IC.vb_0    = quat2mat(body_IC.qb_0)*pathDB.maneuver.az.(char(['wb_' num2str(wb_nr)])).body_kin.uvw(1,:)';
                        body_IC.wb_0    = pathDB.maneuver.az.(char(['wb_' num2str(wb_nr)])).body_kin.wb(1,:)';
                        body_IC.F_0     = DynSim_avg.(char(['seq_' num2str(seq_nr)])).sim_data.F_0;
                        body_IC.M_0     = DynSim_avg.(char(['seq_' num2str(seq_nr)])).sim_data.M_0;
%                         body_IC.F_0     = [0; 0; 0];
%                         body_IC.M_0     = [0; 0; 0];                          
                        
                        [kine_man, sim_data_man] = Dynamic_simulation_man( settings, pathDB, body_IC, 'az', char(az_names(wb_nr)), seq_nr, 1 );
                        
                        DynSim_man.(char(man_wb_names(j))).kine         = kine_man;
                        DynSim_man.(char(man_wb_names(j))).sim_data     = sim_data_man;
                        
                    end
                    
                    if strncmpi(man_wb_names(j),'wx',2) == 1
                        
                        t_string = char(man_wb_names(j));
                        wb_nr = str2num(t_string(7:end));
                        
                        seq_nr          = pathDB.maneuver.wx.(char(['wb_' num2str(wb_nr)])).seq_nr;
                        body_IC.xyz_0   = [0; 0; 0];
                        body_IC.qb_0    = pathDB.maneuver.wx.(char(['wb_' num2str(wb_nr)])).body_kin.qb(1,:)';
                        body_IC.vb_0    = quat2mat(body_IC.qb_0)*pathDB.maneuver.wx.(char(['wb_' num2str(wb_nr)])).body_kin.uvw(1,:)';
                        body_IC.wb_0    = pathDB.maneuver.wx.(char(['wb_' num2str(wb_nr)])).body_kin.wb(1,:)';
                        body_IC.F_0     = DynSim_avg.(char(['seq_' num2str(seq_nr)])).sim_data.F_0;
                        body_IC.M_0     = DynSim_avg.(char(['seq_' num2str(seq_nr)])).sim_data.M_0;
%                         body_IC.F_0     = [0; 0; 0];
%                         body_IC.M_0     = [0; 0; 0];                          
                        
                        [kine_man, sim_data_man] = Dynamic_simulation_man( settings, pathDB, body_IC, 'wx', char(wx_names(wb_nr)), seq_nr, 1 );
                        
                        DynSim_man.(char(man_wb_names(j))).kine         = kine_man;
                        DynSim_man.(char(man_wb_names(j))).sim_data     = sim_data_man;
                        
                    end
                    
                    if strncmpi(man_wb_names(j),'wy',2) == 1
                        
                        t_string = char(man_wb_names(j));
                        wb_nr = str2num(t_string(7:end));
                        
                        seq_nr          = pathDB.maneuver.wy.(char(['wb_' num2str(wb_nr)])).seq_nr;
                        body_IC.xyz_0   = [0; 0; 0];
                        body_IC.qb_0    = pathDB.maneuver.wy.(char(['wb_' num2str(wb_nr)])).body_kin.qb(1,:)';
                        body_IC.vb_0    = quat2mat(body_IC.qb_0)*pathDB.maneuver.wy.(char(['wb_' num2str(wb_nr)])).body_kin.uvw(1,:)';
                        body_IC.wb_0    = pathDB.maneuver.wy.(char(['wb_' num2str(wb_nr)])).body_kin.wb(1,:)';
                        body_IC.F_0     = DynSim_avg.(char(['seq_' num2str(seq_nr)])).sim_data.F_0;
                        body_IC.M_0     = DynSim_avg.(char(['seq_' num2str(seq_nr)])).sim_data.M_0;
%                         body_IC.F_0     = [0; 0; 0];
%                         body_IC.M_0     = [0; 0; 0];                          
                        
                        [kine_man, sim_data_man] = Dynamic_simulation_man( settings, pathDB, body_IC, 'wy', char(wy_names(wb_nr)), seq_nr, 1 );
                        
                        DynSim_man.(char(man_wb_names(j))).kine         = kine_man;
                        DynSim_man.(char(man_wb_names(j))).sim_data     = sim_data_man;
                        
                    end
                    
                    if strncmpi(man_wb_names(j),'wz',2) == 1
                        
                        t_string = char(man_wb_names(j));
                        wb_nr = str2num(t_string(7:end));
                        
                        seq_nr          = pathDB.maneuver.wz.(char(['wb_' num2str(wb_nr)])).seq_nr;
                        body_IC.xyz_0   = [0; 0; 0];
                        body_IC.qb_0    = pathDB.maneuver.wz.(char(['wb_' num2str(wb_nr)])).body_kin.qb(1,:)';
                        body_IC.vb_0    = quat2mat(body_IC.qb_0)*pathDB.maneuver.wz.(char(['wb_' num2str(wb_nr)])).body_kin.uvw(1,:)';
                        body_IC.wb_0    = pathDB.maneuver.wz.(char(['wb_' num2str(wb_nr)])).body_kin.wb(1,:)';
                        body_IC.F_0     = DynSim_avg.(char(['seq_' num2str(seq_nr)])).sim_data.F_0;
                        body_IC.M_0     = DynSim_avg.(char(['seq_' num2str(seq_nr)])).sim_data.M_0;
%                         body_IC.F_0     = [0; 0; 0];
%                         body_IC.M_0     = [0; 0; 0];                          
                        
                        [kine_man, sim_data_man] = Dynamic_simulation_man( settings, pathDB, body_IC, 'wz', char(wz_names(wb_nr)), seq_nr, 1 );
                        
                        DynSim_man.(char(man_wb_names(j))).kine         = kine_man;
                        DynSim_man.(char(man_wb_names(j))).sim_data     = sim_data_man;
                        
                    end
                    
                    i = i+1;
                    
                end
                
        end
            
        save(savefile2,'DynSim_avg','DynSim_steady','DynSim_man')       
            
    end       
    
%     % Save the data:
%     
%     save(savefile2,'DynSim_avg','DynSim_man','DynSim_steady')
    
    % Make validation plots:
    
    Dynamic_Model_Validation_plots( settings, pathDB, DynSim_avg, DynSim_man, DynSim_steady )

    
    % Save the data:
    
    save(savefile1,'DynSim_avg','DynSim_man','DynSim_steady')

end

