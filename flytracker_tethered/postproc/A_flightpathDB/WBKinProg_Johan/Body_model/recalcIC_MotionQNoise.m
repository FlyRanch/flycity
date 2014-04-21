function [Y,tnew] = recalcIC_MotionQNoise(pp,NN,idx,PAR,mode)
% Y = recalcIC_Joint(pp,NN,idx,PAR,mode)
% This function predicts the new joint angles of the fly by comparing them
% with a database of stored joint angles.

persistent p_ %W mu


%Here, I will call the function and load the previous solutions.  I will
%store it so that I can use them with the prediction function called later
switch mode
    case 'load'
        for j = 1:length(idx)
            load([PAR.solutionpath 'fly_' PAR.stub '/fly' num2str(idx(j)) '.mat']);
            numtot = PAR.statedim;
            for k = PAR.flynum
                if k == 1
                    beg = 1;
                else
                    beg = sum(numtot(1:k-1))+1;
                end
                p_(:,j,k) = xh(beg:beg+PAR.statedim-1);
                
            end
        end
        return
        
    case 'calc'
        if isempty(p_)
            error('You have not loaded the previous solutions into memory!')
        end
        Pfull = [p_ pp];   
        [Y,tnew] = predictMotionQNoise(Pfull,NN,PAR);
       
end