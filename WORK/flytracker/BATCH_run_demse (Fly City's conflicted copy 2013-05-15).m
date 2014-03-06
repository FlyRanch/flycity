list = dir('20*');
% list(1:2) = [];

for i = 40:44
    if list(i).isdir == 1
        cd(list(i).name)
        pause(1)
        demse_MODAL_flo_MoveModFlyTracked_AutoEndNkine_flycity
%         demse_MODAL_flo_MoveModFlyTracked_AutoEndNkine_matts
%         demse_MODAL_flo_MoveModFlyTracked_AutoEndNkine_florians
        
        cd ..
    end
end

