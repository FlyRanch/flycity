addpath('/home/florian/Dropbox/WORK/flytracker/flytracker');
dir1=dir
for a=3:length(dir1)
    if dir1(a).isdir==1
        cd(dir1(a).name)
        
        dir2 = dir
        for b=3:length(dir2)
            if dir2(b).isdir==1
                cd(dir2(b).name)
                
                dir_fly = dir('fly*')
                
                if isempty(dir_fly)==1
                    
                elseif dir_fly(1).isdir==0
                    
                    cd ..
                    if exist('Projections')==7
                        cd(dir2(b).name)
                    else
                        cd(dir2(b).name)
                        paste_image_FTMmod_Photronnames_frameloop
                    end
                    
                elseif dir_fly(1).isdir==1
                    
                    cd(dir_fly(1).name)
                    dir_fly2 = dir('fly*')
                    
                    if isempty(dir_fly)==0
                        cd ..
                        if exist('Projections')==7
                            cd(dir_fly(1).name)
                        else
                            cd(dir_fly(1).name)
                            paste_image_FTMmod_Photronnames_frameloop
                        end
                    end
                    cd ..
                end
                cd ..
            end
        end
        cd ..
    end
end
