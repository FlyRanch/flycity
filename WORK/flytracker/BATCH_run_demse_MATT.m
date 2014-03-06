% set working dir to /media/""/solutions/

list = dir;
list(1:2) = [];

for i = 20:22
    if list(i).isdir == 1
        cd(list(i).name)
        pause(1)
        
%         demse_MODAL_flo_MoveModFlyTracked_AutoEndNkine_flycity
          demse_MODAL_flo_MoveModFlyUntracked_AutoEndNkine_MATT
%         demse_MODAL_flo_MoveModFlyTracked_AutoEndNkine_MATT
%         demse_MODAL_flo_MoveModFlyTracked_AutoEndNkine
        
        cd ..
    end
end

    