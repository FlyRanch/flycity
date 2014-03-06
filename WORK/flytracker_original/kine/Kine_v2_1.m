% KINE_v2_1
%
% This program is a tool for extracting 3D coordinates from video data.
% Open the graphical user interface by typing "Kine_v2_0" at the Matlab
% prompt.  It should be able to accomodate a range of camera configurations
% and resolutions, and it should function across computer platforms and
% screens.  The digitized points are grouped into objects of different
% types (points, segments, arcs, models) so that multiple digitization
% techniques may be applied simultaneously.
%
% NOTE: The GUIDE was NOT used to make this GUI, so the program is called
%
% Created by Gwyneth Card, February 9 2004
% Based on a previous version of Kine by Stephen Fry
% With 3D calibration and regression functions written by Will Dickson and 
% Christoph Reinschmidt
%
% 2004.06.28 Version 2_1 includes scripts and debugging by Qing Liu

%---------------------------------------------------------------------
function Kine_v2_1
clear
clc

warning off MATLAB:divideByZero
%---------------------------------------------------------------------
% Look for setup folder in current directory, if not there then load
% hard-wired setup files
if exist('setup') == 7
    set_list = dir([cd,filesep,'setup',filesep,'*.mat']);
    for i = 1:size(set_list,1)
        setups{i} = set_list(i).name;
    end
else
    setups = {'default.mat',...
        'twocam.mat',...
        'hummingbird.mat',...
        'mike.mat','takeoff.mat'};
end

[setup_num, v] = listdlg('PromptString','Choose a Setup file',...
    'SelectionMode','single',...
    'ListString',setups);

setup_file = setups{setup_num};
          
setup_path = [cd,filesep,'setup',filesep];

found_file = exist([setup_path,setup_file]);
if found_file == 2
    data.setup = load([setup_path,setup_file]);
else
    msg = 'No default setup file detected, please create a setup called "default" using the Setup choice in the Options menu';
    dlgname = 'NO DEFAULT SETUP';
    warndlg(msg,dlgname)
    data.setup = [];
end

% Check to see that the default setup file has all the correct fields
setup_wrong = 0;
if isfield(data.setup,'name') == 0
    setup_wrong = setup_wrong + 1;
    data.setup.name = 'hard-wired';
end
if isfield(data.setup,'cam_num') == 0
    setup_wrong = setup_wrong + 1;
    data.setup.cam_num = 3;
end
if isfield(data.setup,'cal_type') == 0
    setup_wrong = setup_wrong + 1;
    data.setup.cal_type = 'DLT';
end
if isfield(data.setup,'cal_mfile') == 0
    setup_wrong = setup_wrong + 1;
    data.setup.cal_mfile = 'get_calibration_points_g';
end
if isfield(data.setup,'image_filter') == 0
    setup_wrong = setup_wrong + 1;
    data.setup.image_filter = '*.bmp';
end
if isfield(data.setup,'mfile_path') == 0
    setup_wrong = setup_wrong + 1;
    data.setup.mfile_path = cd;
elseif exist(data.setup.mfile_path) == 0
    data.setup.mfile_path = uigetdir('','Locate the Kine_v2_0 folder');
end
if isfield(data.setup,'data_path') == 0
    setup_wrong = setup_wrong + 1;
    data.setup.data_path = cd;
elseif exist(data.setup.data_path) == 0
    data.setup.data_path = uigetdir('','Locate the data folder');
end
if isfield(data.setup,'keyboard_command_file') == 0
    setup_wrong = setup_wrong + 1;
    data.setup.keyboard_command_file = 'keyboard_commands';
end

if setup_wrong ~= 0
    msg = ['Hard-wired setup values were used for ',num2str(setup_wrong),' parameter(s).'];
    dlgname = 'HARD-WIRED SETUP VALUES USED';
    warndlg(msg,dlgname)
end

cd(data.setup.mfile_path)

%---------------------------------------------------------------------
% Add path for internal folders with mfiles to be accessed all the time
addpath([data.setup.mfile_path,filesep,'calibration_programs'])
addpath([data.setup.mfile_path,filesep,'kine_math'])

%---------------------------------------------------------------------
% Create save parameters
data.save.pathname = '';
data.save.filename = '';
data.save.timestamp = 'Not Saved';
data.save.autosave = 'off';

%---------------------------------------------------------------------
% Set default advance setup
data.advance = 2;

%---------------------------------------------------------------------
% Create the controller figure window
controller = figure(...
    'Tag','controller',...
    'Units','normalized',...
    'Position',[0 0.05 0.3 .91],...
    'Name','KINE CONTROLLER',...
    'NumberTitle','off',...
    'MenuBar','none',...
    'KeyPressFcn',data.setup.keyboard_command_file,...
    'CloseRequestFcn','controller_closereq'...
    );

%---------------------------------------------------------------------
% Create color palate (could make this part of setup later)

data.colors.red     = [ 1.000 0.000 0.000 ];
data.colors.orange  = [ 1.000 0.502 0.000 ];
data.colors.yellow  = [ 1.000 1.000 0.000 ];
data.colors.green   = [ 0.000 0.502 0.000 ];
data.colors.cyan    = [ 0.000 1.000 1.000 ];
data.colors.blue    = [ 0.000 0.000 1.000 ];
data.colors.purple  = [ 0.502 0.000 0.502 ];
data.colors.magenta = [ 1.000 0.000 1.000 ];
data.colors.lemon   = [ 0.910 0.950 0.650 ];     
data.colors.lime    = [ 0.000 1.000 0.000 ];
data.colors.aqua    = [ 0.000 0.502 0.502 ];
data.colors.sky     = [ 0.000 0.502 1.000 ];
data.colors.pink    = [ 1.000 0.502 1.000 ];
data.colors.black   = [ 0.000 0.000 0.000 ]; 
data.colors.white   = [ 1.000 1.000 1.000 ];
data.colors.gray    = [ 0.850 0.850 0.850 ];       

background = get(controller,'Color');   % This should get the default background color
data.colors.background = background;    % which may be different on different systems
          
button_color = data.colors.lemon;
obj_color = data.colors.lemon;
data.colors.button_color = button_color;
data.colors.obj_color = obj_color;
%---------------------------------------------------------------------
% Get list of object types available from Kine_v2_0\obj_types folder
data.types = get_types(data.setup.mfile_path);
if isempty(data.types)
    warndlg('Object Types folder cannot be found on the specified path','No Object Types folder')
    delete(gcf)
    return
end

%---------------------------------------------------------------------
% Create the controller figure menu
controller_menu = uimenu(...
    'Label','Options',...
    'Tag','controller_menu'...
    );

data_menu = uimenu(...
    'Label','Data',...
    'Tag','data_menu'...
    );

opt_menu_setup = uimenu(controller_menu,...
    'Label','Setup...',...
    'Tag','opt_menu_setup'...
    );

opt_menu_objects = uimenu(controller_menu,...
    'Label','Objects...',...
    'Tag','opt_menu_objects',...
    'Callback','open_object_dialog'...
    );

opt_menu_advance = uimenu(controller_menu,...
    'Label','Auto advance...',...
    'Tag','opt_menu_advance',...
    'Callback','choose_advance_option'...
    );

opt_menu_bcviews = uimenu(controller_menu,...
    'Label','Show body-centered views...',...
    'Tag','opt_menu_bcviews',...
    'Callback','show_bcviews'...
    );

opt_menu_saveas = uimenu(controller_menu,...
    'Label','Save As...',...
    'Separator','on',...
    'Tag','opt_menu_saveas',...
    'Callback','save_as'...
    );

opt_menu_view = uimenu(controller_menu,...
    'Label','View Data',...
    'Tag','opt_menu_view'...
    );

data_menu_spline = uimenu(data_menu,...
    'Label','Spline Data...',...
    'Tag','data_menu_spline',...
    'Callback','spline_gui'...
    );

data_menu_unspline = uimenu(data_menu,...
    'Label','UNspline Data...',...
    'Tag','data_menu_unspline',...
    'Callback','unspline_gui'...
    );

data_menu_filter = uimenu(data_menu,...
    'Label','Filter Data...',...
    'Tag','data_menu_filter',...
    'Callback','filter_gui'...
    );

data_menu_calc_updateBCangles = uimenu(data_menu,...
    'Label','Update Body Centered Angles',...
    'Tag','data_menu_calc_updateBCangles',...
    'Callback','update_body_centered_angles'...
    );
    
data_menu_export_rw = uimenu(data_menu,...
    'Label','Export Real World Angles...',...
    'Tag','data_menu_export_rw',...
    'Callback','export_real_world_angles'...
    );

data_menu_export_ha = uimenu(data_menu,...
    'Label','Export Horizontal Angles...',...
    'Tag','data_menu_export_ha',...
    'Callback','export_horizontal_angles'...
    );

data_menu_export_bc = uimenu(data_menu,...
    'Label','Export Body Centered Angles...',...
    'Tag','data_menu_export_bc',...
    'Callback','export_body_centered_angles'...
    );
          
                                            
%---------------------------------------------------------------------
% Create the frames
load_frame = uicontrol(...
    'Style','frame',...
    'Tag','load_frame',...
    'Units','normalized',...
    'Position',[0.05 0.745 0.9 0.235],...
    'BackgroundColor',background...
    );

get_frame = uicontrol(...
    'Style','frame',...
    'Tag','get_frame',...
    'Units','normalized',...
    'Position',[0.05 0.225 0.9 0.495],...
    'BackgroundColor',background,...
    'Visible','off'...
    );

vert_divide = uicontrol(...
    'Style','frame',...
    'Tag','vert_divide',...
    'Units','normalized',...
    'Position',[.342 .26 .002 .36],...
    'BackgroundColor',background,...
    'Visible','off'...
    );

obj_frame = uicontrol(...
    'Style','frame',...
    'Tag','obj_frame',...
    'Units','normalized',...
    'Position',[0.05 0.02 0.9 0.18],...
    'BackgroundColor',background,...
    'Visible','off'...
    );

load_header = uicontrol(...
    'Style','text',...
    'Tag','load_header',...
    'Units','normalized',...
    'Position',[.13 0.965 0.43 0.022],...
    'BackgroundColor',background,...
    'FontName','Helvetica',...
    'FontUnits','points',...
    'FontSize',12,...
    'FontWeight','bold',...
    'String','Load Workspace'...
    );

get_header = uicontrol(...
    'Style','text',...
    'Tag','get_header',...
    'Units','normalized',...
    'Position',[.13 0.706 0.43 0.022],...
    'BackgroundColor',background,...
    'FontName','Helvetica',...
    'FontUnits','points',...
    'FontSize',12,...
    'FontWeight','bold',...
    'String','Get Kinematics',...
    'Visible','off'...
    );

obj_header = uicontrol(...
    'Style','text',...
    'Tag','obj_header',...
    'Units','normalized',...
    'Position',[.13 0.19 0.43 0.023],...
    'BackgroundColor',background,...
    'FontName','Helvetica',...
    'FontUnits','points',...
    'FontSize',12,...
    'FontWeight','bold',...
    'String','Object Controls',...
    'Visible','off'...
    );

%---------------------------------------------------------------------
% Create the buttons & text
new_text = uicontrol(...
    'Style','text',...
    'Tag','new_text',...
    'Units','normalized',...
    'Position',[0.15 0.92 0.2 0.028],...
    'BackgroundColor',background,...
    'FontName','Helvetica',...
    'FontUnits','points',...
    'FontSize',11,...
    'FontWeight','bold',...
    'String','NEW:'...
    );
    
saved_text = uicontrol(...
    'Style','text',...
    'Tag','saved_text',...
    'Units','normalized',...
    'Position',[0.63 0.92 0.2 0.028],...
    'BackgroundColor',background,...
    'FontName','Helvetica',...
    'FontUnits','points',...
    'FontSize',11,...
    'FontWeight','bold',...
    'String','SAVED:'...
    );    

load_images_button = uicontrol(...
    'Style','pushbutton',...
    'Tag','load_images_button',...
    'Units','normalized',...
    'Position',[0.1 0.89 0.33 0.03],...
    'BackgroundColor',button_color,...
    'FontName','Helvetica',...
    'FontUnits','points',...
    'FontSize',10,...
    'FontWeight','normal',...
    'String','Load Images',...
    'ToolTipString','Start by loading an image sequence to digitize',...
    'Callback','load_images'...
    );

calibrate_button = uicontrol(...
    'Style','pushbutton',...
    'Tag','calibrate_button',...
    'Units','normalized',...
    'Position',[0.1 0.85 0.33 0.03],...
    'BackgroundColor',button_color,...
    'FontName','Helvetica',...
    'FontUnits','points',...
    'FontSize',10,...
    'FontWeight','normal',...
    'String','Calibrate',...
    'ToolTipString','Calibrate the images by loading a file or running a calibration script',...
    'Callback','load_calibration'...
    );

load_data_button = uicontrol(...
    'Style','pushbutton',...
    'Tag','load_data_button',...
    'Units','normalized',...
    'Position',[0.55 0.85 0.33 0.07],...
    'BackgroundColor',button_color,...
    'FontName','Helvetica',...
    'FontUnits','points',...
    'FontSize',10,...
    'FontWeight','normal',...
    'String','Load Data',...
    'ToolTipString','Loads the images, calibration and data for an already analyzed sequence',...
    'Callback','load_data'...
    );

image_text = uicontrol(...
    'Style','text',...
    'Tag','image_text',...
    'Units','normalized',...
    'Position',[0.12 0.80 0.75 0.028],...
    'BackgroundColor',background,...
    'FontName','Helvetica',...
    'FontUnits','points',...
    'FontSize',8,...
    'FontWeight','normal',...
    'FontAngle','italic',...
    'String','< No Image Sequence Loaded >'...
    );

cal_text = uicontrol(...
    'Style','text',...
    'Tag','cal_text',...
    'Units','normalized',...
    'Position',[0.12 0.775 0.75 0.028],...
    'BackgroundColor',background,...
    'FontName','Helvetica',...
    'FontUnits','points',...
    'FontSize',8,...
    'FontWeight','normal',...
    'FontAngle','italic',...
    'String','< No Calibration Loaded >'...
    );

data_text = uicontrol(...
    'Style','text',...
    'Tag','data_text',...
    'Units','normalized',...
    'Position',[0.12 0.75 0.75 0.028],...
    'BackgroundColor',background,...
    'FontName','Helvetica',...
    'FontUnits','points',...
    'FontSize',8,...
    'FontWeight','normal',...
    'FontAngle','italic',...
    'String','< Data Not Saved >'...
    );

frame_slider = uicontrol(...
    'Style','slider',...
    'Tag','frame_slider',...
    'Units','normalized',...
    'Position',[0.1 0.65 0.6 0.033],...
    'BackgroundColor',data.colors.gray,...
    'Interruptible','off',...
    'Visible','off',...
    'Callback','update_images'...
    );

frame_box = uicontrol(...
    'Style','edit',...
    'Tag','frame_box',...
    'Units','normalized',...
    'Position',[0.72 0.65 0.17 0.033],...
    'BackgroundColor','white',...
    'FontName','Helvetica',...
    'FontUnits','points',...
    'FontSize',10,...
    'FontWeight','normal',...
    'String','0',...
    'Interruptible','off',...
    'Visible','off',...
    'Callback','update_images'...
    );

frame_text = uicontrol(...
    'Style','text',...
    'Tag','frame_text',...
    'Units','normalized',...
    'Position',[0.72 0.695 0.17 0.015],...
    'BackgroundColor',background,...
    'FontName','Helvetica',...
    'FontUnits','points',...
    'FontSize',10,...
    'FontWeight','normal',...
    'String','FRAME',...
    'Interruptible','off',...
    'Visible','off'...
    );

pointer_text = uicontrol(...
    'Style','text',...
    'Tag','pointer_text',...
    'Units','normalized',...
    'Position',[0.367 0.442 0.23 0.0175],...
    'BackgroundColor',background,...
    'FontName','Helvetica',...
    'FontUnits','points',...
    'FontSize',8,...
    'FontWeight','normal',...
    'HorizontalAlignment','left',...
    'String','Pointer Function:',...
    'Visible','off'...
    );

pointer_menu = uicontrol(...
    'Style','popupmenu',...
    'Tag','pointer_menu',...
    'Units','normalized',...
    'Position',[0.62 0.43 0.3 0.032],...
    'BackgroundColor','white',...
    'FontName','Helvetica',...
    'FontUnits','points',...
    'FontSize',11,...
    'FontWeight','normal',...
    'String',{'Mark Points','Zoom'},...
    'Value',2,...
    'Visible','off',...
    'Callback','toggle_pointer_fcn(''menu_in'')'...
    );

cur_object_text = uicontrol(...
    'Style','text',...
    'Tag','cur_object_text',...
    'Units','normalized',...
    'Position',[0.367 0.402 0.23 0.0175],...
    'BackgroundColor',background,...
    'FontName','Helvetica',...
    'FontUnits','points',...
    'FontSize',8,...
    'FontWeight','normal',...
    'HorizontalAlignment','left',...
    'String','Current Object:',...
    'Visible','off'...
    );

cur_object_menu = uicontrol(...
    'Style','popupmenu',...
    'Tag','cur_object_menu',...
    'Units','normalized',...
    'Position',[0.62 0.39 0.3 0.032],...
    'BackgroundColor','white',...
    'FontName','Helvetica',...
    'FontUnits','points',...
    'FontSize',11,...
    'FontWeight','normal',...
    'String',{'NONE'},...
    'Visible','off',...
    'Callback','set_cur_pt_menu'...
    );

cur_point_text = uicontrol(...
    'Style','text',...
    'Tag','cur_point_text',...
    'Units','normalized',...
    'Position',[0.367 0.362 0.23 0.0175],...
    'BackgroundColor',background,...
    'FontName','Helvetica',...
    'FontUnits','points',...
    'FontSize',8,...
    'FontWeight','normal',...
    'HorizontalAlignment','left',...
    'String','Current Point:',...
    'Visible','off'...
    );

cur_point_menu = uicontrol(...
    'Style','popupmenu',...
    'Tag','cur_point_menu',...
    'Units','normalized',...
    'Position',[0.62 0.35 0.3 0.032],...
    'BackgroundColor','white',...
    'FontName','Helvetica',...
    'FontUnits','points',...
    'FontSize',11,...
    'FontWeight','normal',...
    'String',{'NONE'},...
    'Visible','off',...
    'Callback','cur_pt_menu_callback'...
    );

obj_params_text = uicontrol(...
    'Style','text',...
    'Tag','obj_params_text',...
    'Units','normalized',...
    'Position',[0.095 0.27 0.3 0.0175],...
    'BackgroundColor',background,...
    'FontName','Helvetica',...
    'FontUnits','points',...
    'FontSize',11,...
    'FontWeight','normal',...
    'HorizontalAlignment','left',...
    'String','Object Parameters:',...
    'Visible','off'...
    );

type_text = uicontrol(...
    'Style','text',...
    'Tag','type_text',...
    'Units','normalized',...
    'Position',[0.367 0.33 0.26 0.0175],...
    'BackgroundColor',background,...
    'FontName','Helvetica',...
    'FontUnits','points',...
    'FontSize',8,...
    'FontWeight','normal',...
    'HorizontalAlignment','left',...
    'String','TYPE:',...
    'Visible','off'...
    );

color_text = uicontrol(...
    'Style','text',...
    'Tag','color_text',...
    'Units','normalized',...
    'Position',[0.655 0.33 0.26 0.0175],...
    'BackgroundColor',background,...
    'FontName','Helvetica',...
    'FontUnits','points',...
    'FontSize',8,...
    'FontWeight','normal',...
    'HorizontalAlignment','left',...
    'String','COLOR:',...
    'Visible','off'...
    );

type_list = uicontrol(...
    'Style','listbox',...
    'Tag','type_list',...
    'Units','normalized',...
    'Position',[0.367 0.26 0.26 0.067],...
    'BackgroundColor',data.colors.gray,...
    'FontName','Helvetica',...
    'FontUnits','points',...
    'FontSize',8,...
    'FontWeight','normal',...
    'HorizontalAlignment','left',...
    'String',data.types,...
    'Visible','off',...
    'Callback','type_list_callback'...
    );

color_list = uicontrol(...
    'Style','listbox',...
    'Tag','color_list',...
    'Units','normalized',...
    'Position',[0.655 0.26 0.26 0.067],...
    'BackgroundColor',data.colors.gray,...
    'FontName','Helvetica',...
    'FontUnits','points',...
    'FontSize',8,...
    'FontWeight','normal',...
    'HorizontalAlignment','left',...
    'String',fieldnames(data.colors),...
    'Visible','off',...
    'Callback','color_list_callback'...
    );

vis_check = uicontrol(...
    'Style','checkbox',...
    'Tag','vis_check',...
    'Units','normalized',...
    'Position',[0.80 0.23 0.125 0.018],...
    'BackgroundColor',background,...
    'FontName','Helvetica',...
    'FontUnits','points',...
    'FontSize',8,...
    'FontWeight','normal',...
    'String','Visible',...
    'Visible','off',...
    'Callback','vis_check_callback'...
    );

%--------------------------------------------------------------------------
reset_button = uicontrol(...
    'Style','pushbutton',...
    'Tag','reset_button',...
    'Units','normalized',...
    'Position',[0.095 0.57 0.215 0.035],...
    'BackgroundColor',obj_color,...
    'FontName','Helvetica',...
    'FontUnits','points',...
    'FontSize',8,...
    'FontWeight','normal',...
    'String','Reset Images',...
    'ToolTipString','Click here to un-zoom images',...
    'Visible','off',...
    'Callback','reset_button_callback'...
    );

obj_button = uicontrol(...
    'Style','pushbutton',...
    'Tag','obj_button',...
    'Units','normalized',...
    'Position',[0.095 0.52 0.215 0.035],...
    'BackgroundColor',obj_color,...
    'FontName','Helvetica',...
    'FontUnits','points',...
    'FontSize',8,...
    'FontWeight','normal',...
    'String','Edit Obj...',...
    'ToolTipString','Click here to add or delete image objects to digitize',...
    'Visible','off',...
    'Callback','open_object_dialog'...
    );

skip_button = uicontrol(...
    'Style','pushbutton',...
    'Tag','skip_button',...
    'Units','normalized',...
    'Position',[0.095 0.47 0.215 0.035],...
    'BackgroundColor',obj_color,...
    'FontName','Helvetica',...
    'FontUnits','points',...
    'FontSize',8,...
    'FontWeight','normal',...
    'String','Skip Point',...
    'ToolTipString','Click here to skip to the next point without recording any data',...
    'Visible','off',...
    'Callback','advance'...
    );

save_button = uicontrol(...
    'Style','pushbutton',...
    'Tag','save_button',...
    'Units','normalized',...
    'Position',[0.095 0.33 0.215 0.035],...
    'BackgroundColor',obj_color,...
    'FontName','Helvetica',...
    'FontUnits','points',...
    'FontSize',8,...
    'FontWeight','bold',...
    'String','SAVE DATA',...
    'ToolTipString','Click here to save data',...
    'Interruptible','off',...
    'Visible','off',...
    'Callback','save_data'...
    );

autosave_check = uicontrol(...
    'Style','checkbox',...
    'Tag','autosave_check',...
    'Units','normalized',...
    'Position',[0.095 0.23 0.155 0.018],...
    'BackgroundColor',background,...
    'FontName','Helvetica',...
    'FontUnits','points',...
    'FontSize',8,...
    'FontWeight','normal',...
    'String','Autosave',...
    'Visible','off',...
    'Callback','save_auto'...
    );

%--------------------------------------------------------------------------
% Sliders to nudge digitized point around

coord_handles = [];
coord_name = {'x' 'y' 'z'};

for i = 1:3
    
    x = 0.3 + i*0.15;
    
    coord_text = uicontrol(...
        'Style','text',...
        'Tag',[coord_name{i},'_text'],...
        'Units','normalized',...
        'Position',[x 0.62 0.08 0.02],...
        'BackgroundColor',background,...
        'FontName','Helvetica',...
        'FontUnits','points',...
        'FontSize',11,...
        'FontWeight','normal',...
        'String',coord_name{i},...
        'Interruptible','off',...
        'Visible','off'...
        );
    
    coord_slider = uicontrol(...
        'Style','slider',...
        'Tag',[coord_name{i},'_slider'],...
        'Units','normalized',...
        'Min',0,...
        'Max',1,...
        'SliderStep',[0.001 0.01],...        
        'Position',[x 0.51 0.08 0.11],...
        'BackgroundColor',data.colors.gray,...
        'Interruptible','off',...
        'Visible','off',...
        'Callback',[coord_name{i},'_slider_callback']...
        );
    
    coord_edit = uicontrol(...
        'Style','edit',...
        'Tag',[coord_name{i},'_edit'],...
        'Units','normalized',...
        'Position',[x-0.006 0.48 0.09 0.02],...
        'BackgroundColor','white',...
        'FontName','Helvetica',...
        'FontUnits','points',...
        'FontSize',8,...
        'FontWeight','normal',...
        'String','0',...
        'Interruptible','off',...
        'Visible','off',...
        'Callback',[coord_name{i},'_edit_callback']...
        );
    
    coord_handles = [coord_handles, coord_text, coord_slider, coord_edit];

end
    
%---------------------------------------------------------------------
% Save handles
data.handles = guihandles(controller);   % get all handles in current figure and save into data.handles structure (listed by Tag)
    % Group handles
    data.handles.load_workspace = [load_frame, load_header, load_images_button, calibrate_button, image_text, cal_text];
    data.handles.slider = [frame_slider, frame_box, frame_text,];
    data.handles.get_kinematics = [get_frame, get_header, vert_divide, data.handles.slider, pointer_text, pointer_menu,...
                                   cur_object_text, cur_object_menu, cur_point_text, cur_point_menu,...
                                   skip_button, reset_button, save_button, obj_button, autosave_check,...
                                   vis_check, coord_handles, color_text, type_text, type_list, color_list];%obj_params_text, 
    data.handles.obj_controls = [obj_frame, obj_header];
    
guidata(controller,data);                % save data structure into figure data structure using GUIDATA