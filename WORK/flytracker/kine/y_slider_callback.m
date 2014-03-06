% Function:       y_slider_callback
% 
% Description:    This function is called when the slider for the y axis
%                 is changed by the user and updates the point that is
%                 currently being edited.
%                 
% Arguments:      None.
%             
% Modifies:       data - structure containing all digitized data
% Returns:        None.
% Error handling: None.
%
% Revision history:   07/06/2004  Qing Liu - initial revision

function y_slider_callback

controller = findobj('Tag','controller');
data = guidata(controller);

slider_value = get(data.handles.y_slider,'Value');
set(data.handles.y_edit,'String',slider_value);

% Get current object,  point, and frame from menus
choice = get(data.handles.cur_object_menu,'Value');
objects = get(data.handles.cur_object_menu,'String');
cur_object = objects{choice};                               % Note that while the cur_object is retrieved as a string
                                                            % because it will be used to retrieve a structure field, the
cur_point = get(data.handles.cur_point_menu,'Value');       % cur_point is retrieved as a number because it will be used
                                                            % to index the rows of a matrix - this number should be okay without
frame = str2num(get(data.handles.frame_box,'String'));      % verification because the menu is made by retrieving the points
                                                            % string, so the name of the point will be data.kine.(object).points{cur_point}
coords = data.kine.(cur_object).data.coords(cur_point,:,frame);

coords(2) = slider_value;
    
store(coords)

calc_it
plot_it
draw_it