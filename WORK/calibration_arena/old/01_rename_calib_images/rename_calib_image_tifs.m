% rename calibration image files (*.tif) for checkerboard recognition prog

files = dir('*.tif')

if isempty(files)==0
    for i=1:length(files)
        old_name = files(i).name;
        new_name = ['pos',old_name(end-4),'_cam',old_name(4),'.tif'];
        copyfile(old_name,new_name);
    end
end
        
    