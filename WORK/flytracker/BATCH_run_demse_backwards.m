list = dir;
list(1:2) = [];

for i = 9:14
    j = length(list) -i +1;
    
    if list(j).isdir == 1
        cd(list(j).name)
        pause(1)
%         demse_MODAL_flo_MoveModFlyTracked_AutoEndNkine
%         demse_MODAL_flo_MoveModFlyTracked_AutoEndNkine_flycity
%         demse_MODAL_flo_MoveModFlyTracked_AutoEndNkine_flyanalysis
%         demse_MODAL_flo_MoveModFlyTracked_AutoEndNkine_florians
        demse_MODAL_flo_MoveModFlyTracked_AutoEndNkine_matts
        
        cd ..
    end
end

    