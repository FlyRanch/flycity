clear
clc
close all

alldirs=dir;
% for i=3:size(alldirs)-2
for i=6:size(alldirs)-2
    if alldirs(i).isdir==1
        cd(alldirs(i).name);
        
        if exist('flytracks')==7
            cd('flytracks')
            allfiles=dir('*.mat')

            for n=2:length(allfiles)
                load(allfiles(n).name)
                path(n-1,:) = xh(1:3)';
            end

            hold on
%             plot3(path(2:end-100,1),path(2:end-100,2),path(2:end-100,3),'-.')
            plot3(path(1:end,1),path(1:end,2),path(1:end,3),'.')
            
            % start
            plot3(path(1,1),path(1,2),path(1,3),'og')
            
            % trigger point: end-.5sec
            if exist('fly001838.mat')==2
                load('fly001838.mat')
                plot3(xh(1),xh(2),xh(3),'or')
            end
            
            clear path
            cd ..
        end
        cd ..
    end
end
axis equal
grid on