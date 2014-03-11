Panel_com('stop');
Panel_com('set_pattern_id',1);
Panel_com('send_gain_bias',[10 0 0 0]);
Panel_com('send_gain_bias',[1 0 0 0]);
Panel_com('set_position', [11 1]);
Panel_com('set_position', [1 1]);
Panel_com('all_off');

Panel_com('start');