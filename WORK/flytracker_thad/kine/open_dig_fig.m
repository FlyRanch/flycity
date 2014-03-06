% Function OPEN_DIG_FIG
% 
% CALLING FUNCTION: load_images
% ACTIONS: Creates the Digitizing window where the images are plotted
% PARENT PROGRAM: Kine_v2_0
% LAST MODIFIED: April 2, 2004 by gwyneth

function open_dig_fig

controller = findobj('Tag','controller');
data = guidata(controller);

dig_fig = findobj('Tag','dig_fig');
if isempty(dig_fig) == 0    % if there is already a dig_fig window
    dig_fig_closereq        % then close it
end

% Make the color scheme of the dig fig same as that of the controller
dig_data.colors = data.colors;
axes_color = dig_data.colors.black;
bkgd = dig_data.colors.background;

dig_fig = figure(...
    'Tag','dig_fig',...
    'Units','normalized',...
    'Position',[0.305 0.05 0.69 .925],...
    'Name','Digitizing Window',...
    'NumberTitle','off',...
    'MenuBar','none',...
    'WindowButtonMotionFcn','crosshair_hotspots',...
    'KeyPressFcn',data.setup.keyboard_command_file,...
    'BackingStore','off',...
    'DoubleBuffer','on',...
    'Renderer','painters',...
    'CloseRequestFcn','dig_fig_closereq'...
    );

% Position values:
a = 0.0280797;
b = 0.518018;
c = 0.463768;
d = 0.461261;
e = 0.520833;
f = 0.0279279;

cam1 = axes(...
    'Units','normalized',...
    'Position',[ a b c d ],...
    'Color',axes_color,...
    'Tag','cam1'...
    );

cam2 = axes(...
    'Units','normalized',...
    'Position',[ e b c d ],...
    'Color',axes_color,...
    'Tag','cam2'...
    );

cam3 = axes(...
    'Units','normalized',...
    'Position',[ a f c d ],...
    'Color',axes_color,...
    'Tag','cam3'...
    );
set([cam1,cam2,cam3],'DrawMode','fast')

% add edit boxes to label axes (no callbacks)
w_label = 0.02;
h_label = 0.017;
offset = 0.002;
label_color = dig_data.colors.gray;

cam1_xlabel = uicontrol(...
    'Style','edit',...
    'Units','normalized',...
    'Position',[a+c-w_label b-(h_label+offset) w_label h_label],...
    'BackgroundColor',label_color...
    );

cam1_ylabel = uicontrol(...
    'Style','edit',...
    'Units','normalized',...
    'Position',[a-(w_label+offset) b+d-h_label w_label h_label],...
    'BackgroundColor',label_color...
    );

cam2_xlabel = uicontrol(...
    'Style','edit',...
    'Units','normalized',...
    'Position',[e+c-w_label b-(h_label+offset) w_label h_label],...
    'BackgroundColor',label_color...
    );

cam2_ylabel = uicontrol(...
    'Style','edit',...
    'Units','normalized',...
    'Position',[e-(w_label+offset) b+d-h_label w_label h_label],...
    'BackgroundColor',label_color...
    );

cam3_xlabel = uicontrol(...
    'Style','edit',...
    'Units','normalized',...
    'Position',[a+c-w_label f-(h_label+offset) w_label h_label],...
    'BackgroundColor',label_color...
    );

cam3_ylabel = uicontrol(...
    'Style','edit',...
    'Units','normalized',...
    'Position',[a-(w_label+offset) f+d-h_label w_label h_label],...
    'BackgroundColor',label_color...
    );


coord_plot = axes(...
    'Tag','coord_plot',...
    'Units','normalized',...
    'Position',[0.52 0.29 0.35 0.18],...
    'Color',bkgd,...
    'FontUnits','points',...
    'FontSize',8,...
    'XGrid','on',...
    'YGrid','on',...
    'ColorOrder',[1 0 0; 0 0.85 0; 0 0 1],...
    'LineStyleOrder','-|-|-',...
    'LineWidth',0.5,...
    'ButtonDownFcn','',...
    'NextPlot','replacechildren'...
    );
xlabel('Frame','FontSize',10,'FontWeight','normal')
title('3D coords of cur point:  X (RED), Y (GREEN), or Z (BLUE)','FontSize',8,'FontWeight','normal')

y_s = 0.42;  % y for topmost edit
y_t = 0.013; % text vertical offset from edit
y_b = 0.03;  % vert distance between edits
x_e = 0.882; % x for edits
x_t = 0.926; % x for texts
edit_l = 0.035;
edit_h = 0.02;
txt_l = 0.065;
txt_h = 0.03;
spc = x_t - x_e - edit_l;

length_check = uicontrol(...
    'Style','checkbox',...
    'Tag','length_check',...
    'Units','normalized',...
    'FontUnits','points',...
    'FontSize',8,...
    'Position',[x_e (y_s+1.7*y_b) (edit_l+txt_l+spc) edit_h],...
    'String','Show Object Length',...
    'Value',1,...
    'HorizontalAlignment','center',...
    'BackgroundColor',bkgd,...
    'Callback','plot_it'...
    );

angles_check = uicontrol(...
    'Style','checkbox',...
    'Tag','angles_check',...
    'Units','normalized',...
    'FontUnits','points',...
    'FontSize',8,...
    'Position',[x_e+0.02 (y_s+y_b) (edit_l+txt_l+spc) edit_h],...
    'String','Show Angles',...
    'HorizontalAlignment','center',...
    'BackgroundColor',bkgd,...
    'Callback','plot_it'...
    );

fr_st_edit = uicontrol(...
    'Style','edit',...
    'Tag','fr_st_edit',...
    'Units','normalized',...
    'FontUnits','points',...
    'FontSize',8,...
    'Position',[x_e y_s edit_l edit_h],...
    'String','0',...
    'BackgroundColor','w',...
    'Callback','set_plot_vals'...
    );

fr_st_text = uicontrol(...
    'Style','text',...
    'Tag','fr_st_text',...
    'Units','normalized',...
    'FontUnits','points',...
    'FontSize',8,...
    'Position',[x_t y_s-y_t txt_l txt_h],...
    'String','Start Frame',...
    'HorizontalAlignment','left',...
    'BackgroundColor',bkgd...
    );

fr_end_edit = uicontrol(...
    'Style','edit',...
    'Tag','fr_end_edit',...
    'Units','normalized',...
    'FontUnits','points',...
    'FontSize',8,...
    'Position',[x_e y_s-y_b edit_l edit_h],...
    'String','0',...
    'BackgroundColor','w',...
    'Callback','set_plot_vals'...
    );

fr_end_text = uicontrol(...
    'Style','text',...
    'Tag','fr_end_text',...
    'Units','normalized',...
    'FontUnits','points',...
    'FontSize',8,...
    'Position',[x_t y_s-y_t-y_b txt_l txt_h],...
    'String','End Frame',...
    'HorizontalAlignment','left',...
    'BackgroundColor',bkgd...
    );

y_min_edit = uicontrol(...
    'Style','edit',...
    'Tag','y_min_edit',...
    'Units','normalized',...
    'FontUnits','points',...
    'FontSize',8,...
    'Position',[x_e y_s-2*y_b edit_l edit_h],...
    'String','0',...
    'BackgroundColor','w',...
    'Callback','set_plot_vals'...
    );

y_min_text = uicontrol(...
    'Style','text',...
    'Tag','y_min_text',...
    'Units','normalized',...
    'FontUnits','points',...
    'FontSize',8,...
    'Position',[x_t y_s-y_t-2*y_b txt_l txt_h],...
    'String','y-axis min',...
    'HorizontalAlignment','left',...
    'BackgroundColor',bkgd...
    );

y_max_edit = uicontrol(...
    'Style','edit',...
    'Tag','y_max_edit',...
    'Units','normalized',...
    'FontUnits','points',...
    'FontSize',8,...
    'Position',[x_e y_s-3*y_b edit_l edit_h],...
    'String','0',...
    'BackgroundColor','w',...
    'Callback','set_plot_vals'...
    );

y_max_text = uicontrol(...
    'Style','text',...
    'Tag','y_max_text',...
    'Units','normalized',...
    'FontUnits','points',...
    'FontSize',8,...
    'Position',[x_t y_s-y_t-3*y_b txt_l txt_h],...
    'String','y-axis max',...
    'HorizontalAlignment','left',...
    'BackgroundColor',bkgd...
    );

norm_check = uicontrol(...
    'Style','checkbox',...
    'Tag','norm_check',...
    'Units','normalized',...
    'FontUnits','points',...
    'FontSize',8,...
    'Position',[x_e+0.02 (y_s-y_t-3.5*y_b) (edit_l+txt_l+spc) edit_h],...
    'String','Normalized',...
    'HorizontalAlignment','center',...
    'BackgroundColor',bkgd,...
    'Callback','plot_it'...
    );

reset_plot = uicontrol(...
    'Style','pushbutton',...
    'Tag','reset_plot',...
    'Units','normalized',...
    'FontUnits','points',...
    'FontSize',8,...
    'Position',[x_e+0.02 (y_s-y_t-4.3*y_b) 0.08 edit_h],...
    'String','Reset Plot',...
    'BackgroundColor',[.85 .85 .85],...
    'Callback','reset_plot'...
    );

obj_length = axes(...
    'Tag','obj_length',...
    'Units','normalized',...
    'Position',[0.85 0.36 0.01 0.10],...
    'Color',bkgd,...
    'FontUnits','points',...
    'FontSize',8,...
    'Color',axes_color,...
    'XColor',bkgd,...
    'NextPlot','replacechildren',...
    'DrawMode','fast'...
    );

%--------------------------------------------------------------------------
% body centered plots

body_plot_side = axes(...
    'Tag','body_plot_side',...
    'Units','normalized',...
    'Position',[0.52 0.14 0.10 0.10],...
    'Color',bkgd,...
    'FontUnits','points',...
    'FontSize',8,...
    'Color',axes_color,...
    'NextPlot','add',...
    'Visible','off'...
    );

body_plot_rear = axes(...
    'Tag','body_plot_rear',...
    'Units','normalized',...
    'Position',[0.66 0.14 0.10 0.10],...
    'Color',bkgd,...
    'FontUnits','points',...
    'FontSize',8,...
    'Color',axes_color,...
    'NextPlot','add',...
    'Visible','off'...
    );

body_plot_top = axes(...
    'Tag','body_plot_top',...
    'Units','normalized',...
    'Position',[0.80 0.14 0.10 0.10],...
    'Color',bkgd,...
    'FontUnits','points',...
    'FontSize',8,...
    'Color',axes_color,...
    'NextPlot','add',...
    'Visible','off'...
    );

axis1_hist = axes(...
    'Tag','axis1_hist',...
    'Units','normalized',...
    'Position',[0.94 0.14 0.01 0.10],...
    'Color',bkgd,...
    'FontUnits','points',...
    'FontSize',8,...
    'Color',axes_color,...
    'XColor',bkgd,...
    'NextPlot','replacechildren',...
    'Visible','off'...
    );

axis2_hist = axes(...
    'Tag','axis2_hist',...
    'Units','normalized',...
    'Position',[0.97 0.14 0.01 0.10],...
    'Color',bkgd,...
    'FontUnits','points',...
    'FontSize',8,...
    'Color',axes_color,...
    'XColor',bkgd,...
    'NextPlot','replacechildren',...
    'Visible','off'...
    );

side_view_text = uicontrol(...
    'Style','text',...
    'Tag','side_view_text',...
    'Units','normalized',...
    'FontUnits','points',...
    'FontSize',8,...
    'Position',[0.52 0.11 0.10 0.01],...
    'String','SIDE',...
    'HorizontalAlignment','center',...
    'BackgroundColor',bkgd,...
    'Visible','off'...
    );

rear_view_text = uicontrol(...
    'Style','text',...
    'Tag','rear_view_text',...
    'Units','normalized',...
    'FontUnits','points',...
    'FontSize',8,...
    'Position',[0.66 0.11 0.1 0.01],...
    'String','REAR',...
    'HorizontalAlignment','center',...
    'BackgroundColor',bkgd,...
    'Visible','off'...
    );

top_view_text = uicontrol(...
    'Style','text',...
    'Tag','top_view_text',...
    'Units','normalized',...
    'FontUnits','points',...
    'FontSize',8,...
    'Position',[0.80 0.11 0.10 0.01],...
    'String','TOP',...
    'HorizontalAlignment','center',...
    'BackgroundColor',bkgd,...
    'Visible','off'...
    );

ax1_text = uicontrol(...
    'Style','text',...
    'Tag','ax1_text',...
    'Units','normalized',...
    'FontUnits','points',...
    'FontSize',8,...
    'Position',[0.935 0.11 0.02 0.01],...
    'String','AX1',...
    'HorizontalAlignment','center',...
    'BackgroundColor',bkgd,...
    'Visible','off'...
    );

ax2_text = uicontrol(...
    'Style','text',...
    'Tag','ax2_text',...
    'Units','normalized',...
    'FontUnits','points',...
    'FontSize',8,...
    'Position',[0.965 0.11 0.02 0.01],...
    'String','AX2',...
    'HorizontalAlignment','center',...
    'BackgroundColor',bkgd,...
    'Visible','off'...
    );
%--------------------------------------------------------------------------

dig_data.handles = guihandles(dig_fig);
dig_data.temp = [];

% Collect all bc handles:
dig_data.handles.bc_view = [body_plot_side, body_plot_rear, body_plot_top,...
    axis1_hist, axis2_hist, side_view_text, rear_view_text, top_view_text,...
    ax1_text, ax2_text];

guidata(dig_fig,dig_data)