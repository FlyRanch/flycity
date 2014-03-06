%% Make MoveModel from successful run of Demse
% working directory /media/"disk"/solutions/

clc
clear all
% 
% home = '/media/Matt01/solutions/';
% 
% list = dir('20*');

% for i = 1:length(list)
%     
%     if list(i).isdir == 1 % if (i) is a directory
%         cd(list(i).name) % then open it
%         
%         if exist('flytracks') == 7 % if flytracks exists
%             cd flytracks % open it
%             list2 = dir('fly*.mat')
%         
%               for j = 1:length(list2)
%               
%                load(list2(j).name);
%                MoveModel.solnQ(:,i)=xh;
%                                        
%               end
%               
%         cd .. % jump out of flytracks
%         save([cd,'_successful_flytrack_MoveModel.mat'],'MoveModel')
%         end
%         cd(home) % jump out of directory
%     end
% end
%   


% list = dir('20*');
% for i = [1:8 18 21:23 20:32]  
%     
%         if list(i).isdir == 1 % if (i) is a directory
%         cd(list(i).name) % then open it
clear all

        cd flytracks

        list2 = dir('fly*');
 
            for j = 1:length(list2)
                load(list2(j).name);
                MoveModel.solnQ(:,j)=xh;
            end
            
            cd ..
            save([cd,'_successful_flytrack_MoveModel.mat'],'MoveModel')
            cd ..
%         end
% end


