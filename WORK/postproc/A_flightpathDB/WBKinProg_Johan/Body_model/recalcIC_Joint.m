function [Y,tnew] = recalcIC_Joint(pp,idx,PAR,mode)
% Y = recalcIC_Joint(pp,len,PAR,NN,fishnum)
% This function predicts the new joint angles of the fly by comparing them
% with a database of stored joint angles.

persistent p_ %W mu


%Here, I will call the function and load the previous solutions.  I will
%store it so that I can use them with the prediction function called later
switch mode
    case 'load'
        for j = 1:length(idx)
            %% photron structure FTMmod 20120603
            load([PAR.solutionpath PAR.solutiondirname '/fly' sprintf(['%0',num2str(PAR.digits),'d'], idx(j)) '.mat']);
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
        
%% OLD
%         %[Y,tnew] = predictJointLoc(Pfull);
%         [Y,tnew] = predictMotionQ_NoBodyScale(Pfull);
%         %[Y,tnew] = predictJointQ(Pfull);
%         %[Y,tnew] = predictJointBodyQ(Pfull);

%% FTMmod NO rescaling, auto load
%         [Y,tnew] = predictMotionQ_NoBodyScale(Pfull);
% 
%% FTMmod NO rescaling, auto load
        [Y,tnew] = predictMotionQ_NoBodyScale_autoload(Pfull,PAR);

%% FTMmod NO rescaling, Dm @ 7.5k fps
%         [Y,tnew] = predictMotionQ_NoBodyScale_Dm_7500fps_straight(Pfull);
% 
%% FTMmod add rescaling of wingmovement
%         [Y,tnew] = predictMotionQ_WingsOnly(Pfull);
%         
%% FTMmod add rescaling of wing&body movement
%         [Y,tnew] = predictMotionQ_BodyScale(Pfull,PAR);

        %keyboard
end