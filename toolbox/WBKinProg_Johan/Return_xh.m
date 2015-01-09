function Return_xh(settings, pathDB)
    
    file_name = 'fly';
    dirname = 'flytracks';
    %location_dir = 'D:/Filtered_Flytracks/';
    %location_dir = 'D:/body_smoothing_tests/Filtered_corr_F/';
    %location_dir = 'D:/body_smoothing_tests/Filtered_corr_J/';
    %location_dir = 'D:/body_smoothing_tests/Filtered_no_corr/';
    location_dir = 'D:/body_smoothing_tests/Unfiltered/';
    

    % Program which converts the filtered xh back into the .mat files
    alldirs=dir;
    
%    for j=1:size(pathDB.x,2)
    for j=1:1
        
        if alldirs(j+2).isdir==1
        
        cd(alldirs(j+2).name);
        
        mkdir('filtered_flytracks')
        
                if exist(dirname)==7
                    cd(dirname)

                start = find(isnan(pathDB.x(:,j))==0, 1 );
                stop = find(isnan(pathDB.x(:,j))==0, 1, 'last' );
                
                allfiles=dir;
                
                    for i = 1:(stop-start+1)

                        x = pathDB.x_filt(start-1+i,j);
                        y = pathDB.y_filt(start-1+i,j);
                        z = pathDB.z_filt(start-1+i,j);

                        q1b = pathDB.qb1_filt(start-1+i,j);
                        q2b = pathDB.qb2_filt(start-1+i,j);
                        q3b = pathDB.qb3_filt(start-1+i,j);
                        q4b = pathDB.qb4_filt(start-1+i,j);

                        q1lw = pathDB.qL1_filt1(start-1+i,j);
                        q2lw = pathDB.qL2_filt1(start-1+i,j);
                        q3lw = pathDB.qL3_filt1(start-1+i,j);
                        q4lw = pathDB.qL4_filt1(start-1+i,j);

                        q1rw = pathDB.qR1_filt1(start-1+i,j);
                        q2rw = pathDB.qR2_filt1(start-1+i,j);
                        q3rw = pathDB.qR3_filt1(start-1+i,j);
                        q4rw = pathDB.qR4_filt1(start-1+i,j);

                        xh_new = [x; y; z; q1b; q2b; q3b; q4b; q1lw; q2lw; q3lw; q4lw; q1rw; q2rw; q3rw; q4rw];
                        
                        % Obtain file name and InternalVariables
                        name_now = allfiles(i+3).name;
                        
                        Int_var = [];
                        file_name_now = '?';
                        
                        if exist(name_now)==2
                            load(name_now)

                            Int_var = InternalVariablesDS;

                            file_name_now = name_now;
                            
                            clear name_now
                        end
                        
                        savefile = strcat(location_dir , alldirs(j+2).name, '/filtered_flytracks/', file_name_now);
                        
                        InternalVariablesDS = Int_var;
                        
                        xh = xh_new;
                        
                        save(savefile,'InternalVariablesDS','xh')

                    end
                    
                    cd ..

                end
                
                cd ..

        end
        
        
        
    end

end

