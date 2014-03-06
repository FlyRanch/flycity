function modelparamtest

num_params = 15;
data.colors.button_color = 'w';
data.colors.background = [.8 .8 .8];

cur_obj = 'test_obj';
frame = 1;


horiz_marg = 7;
vert_marg = 1;
horiz_btwn = 3;
vert_btwn = 1;
l_edit = 9;
h_edit = 1.6;
l_slide = 2*l_edit;
h_slide = h_edit;
text_gap = 2;
l_text = 15;
h_text = h_edit;


figure
for i = 1:num_params
    
    if i <= 5
        x = horiz_marg;
        y = vert_marg + (5 - i)*(h_edit + vert_btwn);
        
    elseif i > 5 & i <= 10
        x = horiz_marg + 1*(l_slide + l_edit + text_gap + l_text + horiz_btwn);
        y = vert_marg + (10 - i)*(h_edit + vert_btwn);
        
    elseif i > 10 & i <= 15
        x = horiz_marg + 2*(l_slide + l_edit + text_gap + l_text + horiz_btwn);
        y = vert_marg + (15 - i)*(h_edit + vert_btwn);
        
    end
    
    
    uicontrol(...
        'Style','slider',...
        'Tag','obj_control',...
        'Units','characters',...
        'Position',[x y l_slide h_slide],...
        'BackgroundColor',data.colors.button_color,...
        'FontName','Helvetica',...
        'FontUnits','points',...
        'FontSize',12,...
        'FontWeight','normal',...
        'String','-',...
        'UserData',i,...
        'Callback','param_slider_callback'...
        );

    uicontrol(...
        'Style','edit',...
        'Tag','obj_control',...
        'Units','characters',...
        'Position',[x+l_slide y l_edit h_edit],...
        'BackgroundColor','white',...
        'FontName','Helvetica',...
        'FontUnits','points',...
        'FontSize',12,...
        'FontWeight','normal',...
        'String','string1',...
        'UserData',i,...
        'Callback','param_edit_callback'...
        );

    uicontrol(...
        'Style','text',...
        'Tag','obj_control',...
        'Units','characters',...
        'Position',[x+l_slide+l_edit+text_gap y-0.2 l_text h_text],...
        'BackgroundColor',data.colors.background,...
        'FontName','Helvetica',...
        'FontUnits','points',...
        'FontSize',11,...
        'FontWeight','normal',...
        'HorizontalAlignment','left',...
        'UserData',i,...
        'String','string2'....
        );    
    
end