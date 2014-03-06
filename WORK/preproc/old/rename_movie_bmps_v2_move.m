% rename calibration image files (*.tif) for checkerboard recognition prog

files = dir('*.bmp')

if isempty(files)==0
    for i=1:length(files)
        old_name = files(i).name;
        new_name = ['frame',old_name(end-6:end)];
%         copyfile(old_name,new_name);
        movefile(old_name,new_name);
    end
end
