% make top projection movie for all seqs
clear
clc

skip = 10
cam = 3

dir1=dir
for a=3:length(dir1)
    if dir1(a).isdir==1
        cd(dir1(a).name)
        
        dir2 = dir('flytracks*')
        for b=1:length(dir2)
            if dir2(b).isdir==1
                cd(dir2(b).name)
                
                dir_fly = dir('fly*')
                
                if isempty(dir_fly)==0
                    paste_top_projection_makemovie
                end
                cd ..
            end
        end
        cd ..
    end
end
