% rename calibration image files (*.tif) for checkerboard recognition prog
clear;
clc;

loc = cd;
dirs = dir;
for j=3:length(dirs)
    if dirs(j).isdir==1
        old_dir = dirs(j).name;
        new_dir = ['cam',old_dir(2:4)];
        mkdir(new_dir);
        cd(old_dir)
        
        % info file
        files = dir('*.cih');
        if isempty(files)==0
            for i=1:length(files)
                    old_name = [loc,'/',old_dir,'/',files(i).name];
                    new_name = [loc,'/',new_dir,'/settings.cih'];
%                     copyfile(old_name,new_name);
                    movefile(old_name,new_name);
            end
        end
        
        % movie frames
        files = dir('*.bmp');
        if isempty(files)==0
            for i=1:length(files)
                    old_name = [loc,'/',old_dir,'/',files(i).name];
                    new_name = [loc,'/',new_dir,'/frame',old_name(end-6:end)];
%                     copyfile(old_name,new_name);
                    movefile(old_name,new_name);
            end
        end
        
        cd ..
    end
end
            
