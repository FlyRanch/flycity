% make top flighttrack movie for all seqs
clear
clc

dir1=dir
for a=3:length(dir1)
    if dir1(a).isdir==1
        cd(dir1(a).name)
        
        dir2 = dir('flytracks*')
        for b=1:length(dir2)
            if dir2(b).isdir==1
                    plot_singletrack_v2_angles
            end
        end
        cd ..
    end
end
