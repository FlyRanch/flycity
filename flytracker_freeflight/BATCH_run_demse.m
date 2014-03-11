list = dir('20*');
% list(1:2) = [];

for i = 18:20
    if list(i).isdir == 1
        cd(list(i).name)
        pause(1)
%         demse_MODAL_flo_MoveModFlyTracked_AutoEndNkine_flycity
%         demse_MODAL_flo_MoveModFlyTracked_AutoEndNkine_flyanalysis
        demse_MODAL_flo_MoveModFlyTracked_AutoEndNkine_matts
%         demse_MODAL_flo_MoveModFlyTracked_AutoEndNkine_florians
        close all
        cd ..
    end
end

