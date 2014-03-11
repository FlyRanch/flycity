function read_raw_data(settings,pathDB)

    savefile = 'pathDB1.mat';
    file_name = 'fly';
    dirname = 'flytracks';

    digits = 6;
    start_dir = 3;
    start_file = 4; % exclude ManualFit
    
    if exist('pathDB') == 1
        
    t =     pathDB.t;
    x =     [];
    y =     [];
    z =     [];
    qB1 =   [];
    qB2 =   [];
    qB3 =   [];
    qB4 =   [];
    qL1 =   [];
    qL2 =   [];
    qL3 =   [];
    qL4 =   [];
    qR1 =   [];
    qR2 =   [];
    qR3 =   [];
    qR4 =   [];
    
    end
    
    
    seq = 0;
    
    cd(char(settings.data_loc))
    
    alldirs=dir;
    
    temp_foldernames = [];
    
    for i=start_dir:length(alldirs)
        
        length(alldirs)-i+1
                                
        temp_foldernames = [temp_foldernames; alldirs(i).name ];
        
        if alldirs(i).isdir == 1
            
            cd(alldirs(i).name);
            
            

            if exist(dirname)==7
                
                cd(dirname)

                seq=seq+1;
                settings.seq(seq,1) = str2num(alldirs(i).name(1:8));
                settings.seq(seq,2) = str2num(alldirs(i).name(11:14));
                
                x(1:length(t),seq) =    nan;
                y(1:length(t),seq) =    nan;
                z(1:length(t),seq) =    nan;
                qB1(1:length(t),seq) =  nan;
                qB2(1:length(t),seq) =  nan;
                qB3(1:length(t),seq) =  nan;
                qB4(1:length(t),seq) =  nan;
                qL1(1:length(t),seq) =  nan;
                qL2(1:length(t),seq) =  nan;
                qL3(1:length(t),seq) =  nan;
                qL4(1:length(t),seq) =  nan;
                qR1(1:length(t),seq) =  nan;
                qR2(1:length(t),seq) =  nan;
                qR3(1:length(t),seq) =  nan;
                qR4(1:length(t),seq) =  nan;

                allfiles=dir;
                
                for n=start_file:length(allfiles)
                    
                    name_now = allfiles(n).name;
                    
                    if exist(name_now)==2
                        
                        load(name_now)

                        frame_now = str2num(name_now(length(file_name)+1:length(file_name)+digits));
                        
                        x(frame_now,seq) =    xh(1);
                        y(frame_now,seq) =    xh(2);
                        z(frame_now,seq) =    xh(3);
                        qB1(frame_now,seq) =  xh(4);
                        qB2(frame_now,seq) =  xh(5);
                        qB3(frame_now,seq) =  xh(6);
                        qB4(frame_now,seq) =  xh(7);
                        qL1(frame_now,seq) =  xh(8);
                        qL2(frame_now,seq) =  xh(9);
                        qL3(frame_now,seq) =  xh(10);
                        qL4(frame_now,seq) =  xh(11);
                        qR1(frame_now,seq) =  xh(12);
                        qR2(frame_now,seq) =  xh(13);
                        qR3(frame_now,seq) =  xh(14);
                        qR4(frame_now,seq) =  xh(15);
                        
                        
                    end
                    
                                                  
                end
                cd ..
            end
            cd ..
        end
    end
    
    xyz =   nan(length(t),3,seq);
    qB =    nan(length(t),4,seq);
    qL =    nan(length(t),4,seq);
    qR =    nan(length(t),4,seq);
    
       
    for m = 1:seq
        
        xyz(:,:,m) = [x(:,m) y(:,m) z(:,m)];
        qB(:,:,m) = [qB1(:,m) qB2(:,m) qB3(:,m) qB4(:,m)];
        qL(:,:,m) = [qL1(:,m) qL2(:,m) qL3(:,m) qL4(:,m)];
        qR(:,:,m) = [qR1(:,m) qR2(:,m) qR3(:,m) qR4(:,m)];
        
    end
    
    sequence_names = temp_foldernames;
    nr_of_seq = seq;
   
    
   
    start_stop = nan(seq,2);
    
    for m = 1:seq
        
        start_stop(m,1) = find(isnan(xyz(:,1,m))==0, 1 );
        start_stop(m,2) = find(isnan(xyz(:,1,m))==0, 1, 'last' );
        
    end
    
    save(savefile,'t','xyz','qB','qL','qR','sequence_names','nr_of_seq','start_stop')

end

