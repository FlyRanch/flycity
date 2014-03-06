% Function SAVE_DATA
% 
% CALLING FUNCTION: Callback for save button, called by save_as
% ACTIONS: Saves variables 'data' and 'timestamp' (created here) into the
%          folder designated by data.setup.save_path
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: March 23, 2004 by gwyneth

function save_data

controller = findobj('Tag','controller');
data = guidata(controller);
cd(data.setup.data_path)

while strcmp(data.save.pathname,'')
    %% Photron structure FTMmod20120611
    if length(data.images.path)>40
        save_path = [data.images.path(1:end-36),'solutions/',data.images.path(end-22:end-15),'_',data.images.path(end-5:end-1)];
    else
        save_path = [];
    end
    if exist(save_path)==7
       cd(save_path)
    else
       mkdir(save_path)
       cd(save_path)
    end
    %%
    data.save.pathname = uigetdir('','Locate the folder in which you would like to save the data');
    if data.save.pathname == 0, return, end
end
while strcmp(data.save.filename,'')
    
    %% Photron structure FTMmod20120611
    if length(data.images.path)>23
        save_name = ['kine_',data.images.path(end-22:end-15),'_',data.images.path(end-5:end-1)];
    end
    data.save.filename = inputdlg('Filename: ','Save',1,{save_name});
%     data.save.filename = inputdlg('Filename: ','Save',1);
    data.save.filename = data.save.filename{1};
    if isempty(data.save.filename), return, end
end
    
if exist(data.save.pathname) == 0
    warndlg('Can''t find folder data was previously saved in, try ''Save As...''', 'Save Path Invalid')
    return
end
    
save_name = [data.save.pathname,filesep,data.save.filename];
data.save.timestamp = datestr(now); % current date and time of save

% Save the data
save(save_name, 'data')

% Update the corresponding text
set_text('data_text',['Last save to ''',data.save.filename,''' at ',data.save.timestamp])

guidata(controller,data); % save chosen path and filenames into structure