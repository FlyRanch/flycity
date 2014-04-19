function [Pnew,tnew] = predictMotionQ_NoBodyScale_autoload(Pfull,PAR)

% This function predicts the n+1 state given the previous m solutions
% See section 3.1 of B. Rosenhahn, "Scaled Motion Dynamics for Markerless      
% Motion Capture", CVPR 2007.
% 
% Right now it calculates by brute force enumeration of all possible
% values.  Possibly a better way to do this.
%
% This function predicts the new body transformation using a translation
% vector and body quaternion
%
% the joint angles are predicted using the quaternion representation

persistent data dataQ transidx

% Joint indices
jidx = 8:15;

P = Pfull(jidx,:);
% load data if it is not already in memory

m = size(P,2);

%%
if isempty(data)
    %% FTMmod movement model
%     TSOL=[];
%     save('TSOL.mat','TSOL')
    
    % autoload multiple files from PAR
    data=[];
    dataQ=[];
    transidx=[];
    
    for i=1:length(PAR.MoveModelNameS)
        load([PAR.MoveModelPath,PAR.MoveModelNameS{i}])
        
        if isfield(MoveModel,'solnQ')
            data = [data MoveModel.solnQ(1:3,:)];
            dataQ = [dataQ MoveModel.solnQ(4:15,:)];
            transidx(i) = size(MoveModel.solnQ,2);
        else
            data = [data MoveModel.stateQ(1:3,:)];
            dataQ = [dataQ MoveModel.stateQ(4:15,:)];
            transidx(i) = size(MoveModel.stateQ,2);
        end
    end
    %THL mod 20140218 eliminate all dynamics from body model.
    data = repmat(data(:,1)',length(data),1)';
    dataQ(1,:) = dataQ(1,1);
    dataQ(2,:) = dataQ(2,1);
    dataQ(3,:) = dataQ(3,1);
    dataQ(4,:) = dataQ(4,1);
    %keyboard
    transidx = cumsum(transidx);
    transidx = [m-1 transidx]; % add m-1 to transidx
end
        
%% single file
% if isempty(data)
%     if isfield(MoveModel,'solnQ')
%         data = [MoveModel.solnQ(1:3,:)];
%         dataQ = [MoveModel.solnQ(4:15,:)];
%         transidx = [size(MoveModel.solnQ,2)];
%         transidx = cumsum(transidx);
%     else
%         data = [MoveModel.stateQ(1:3,:)];
%         dataQ = [MoveModel.stateQ(4:15,:)];
%         transidx = [size(MoveModel.stateQ,2)];
%         transidx = cumsum(transidx);
%     end
% end
%% OLD    
% if isempty(data)
%     load('exp101_ManualTrack.mat')
%     load('exp083_ManualTrack.mat')
% %     load('results/exp093_ManualTrack.mat')
% %     keyboard
%     data = [StateQ_sym(1:3,:) StartleStateQ_sym(1:3,:)];
%     dataQ = [StateQ_sym(4:15,:) StartleStateQ_sym(4:15,:)];
%     transidx = [size(StateQ_sym,2) size(StartleStateQ_sym,2)];
%     
% %     data = [StateQ_sym(1:3,:) StartleStateQ_sym(1:3,:) WideSweepStateQ(1:3,:)];
% %     dataQ = [StateQ_sym(4:15,:) StartleStateQ_sym(4:15,:) WideSweepStateQ(4:15,:)];
% %     transidx = [size(StateQ_sym,2) size(StartleStateQ_sym,2) size(WideSweepStateQ,2)];
%     transidx = cumsum(transidx);
% end
%%

%this is the range of joint angle scalings
% scale = [8/12 8/10 1 1.2 4/3];
% scale = 1;
% scale = .95:.01:1.05;
scale = PAR.MoveModelScale;

tt = [1 size(data,2)];

minval = zeros(1,length(scale));
minidx = zeros(1,length(scale));

AllData = cell(length(scale),1);
AllTime = cell(length(scale),1);

for k = 1:length(scale)
    
    scaledtime = tt(1):1/scale(k):tt(2);
    
    %initialize
    scaleddataQ = zeros(size(dataQ,1),length(scaledtime));
    %Scale the body and joint quaternion variables
    scaleddataQ(1:4,:) = slerpMAT(dataQ(1:4,:),scale(k));
    scaleddataQ(5:8,:) = slerpMAT(dataQ(5:8,:),scale(k));
    scaleddataQ(9:12,:) = slerpMAT(dataQ(9:12,:),scale(k));
    
    %scale the body translation data
    scaleddata = zeros(size(data,1),length(scaledtime));
    for i = 1:size(data,1)
        scaleddata(i,:) = interp1(tt(1):tt(2),data(i,:),scaledtime);
    end
    
    AllData{k} = [scaleddata ;scaleddataQ];
    AllTime{k} = scaledtime;
    
    Error = zeros(1,length(scaledtime));
    for i = m:size(scaleddataQ,2)
        Error(i) = sum( sqrt( sum((scaleddataQ(5:12,i-m+1:i) - P).^2,1) ) );
    end
    

    
%% FTMmod: find frame with minimum error EXCLUDING first m-1 (zero) values and last point from each seq (transidx)
%     scaledtransidx = round(transidx*scale(k));
%     clear minval_tmp minidx_tmp
%     for i = 1:length(scaledtransidx)-1
%         [minval_tmp(i),minidx_tmp(i)] = min(Error(scaledtransidx(i)+1:scaledtransidx(i+1)-1));
%         minidx_tmp(i) = minidx_tmp(i) + scaledtransidx(i);
%     end
%     minidx(k) = minidx_tmp(find(minval_tmp==min(minval_tmp)));
%     minval(k) = min(minval_tmp);
% 

% !!! this does not work for scaled data (scale ~= 1) !!!
    clear minval_tmp minidx_tmp
    for i = 1:length(transidx)-1
        [minval_tmp(i),minidx_tmp(i)] = min(Error(transidx(i)+1:min(transidx(i+1)-1,size(Error,2)-1)));
        minidx_tmp(i) = minidx_tmp(i) + transidx(i);
    end
    minidx(k) = min(minidx_tmp(find(minval_tmp==min(minval_tmp))));
    minval(k) = min(minval_tmp);    
% 
%     [minval(k),minidx(k)] = min(Error(m:end-1));
%     %fix the index because we ignore the first m-1 error values because
%     %they are zero by construction.  Also ignore the last value in case it
%     %is the prediction
%     minidx(k) = minidx(k) + (m-1);
end

% scaleidx is the scale with the minimum error
[val,scaleidx] = min(minval);

% This is the n+1 state that is predicted
tsol = minidx(scaleidx);

% load('TSOL.mat')
% TSOL(end+1)=tsol;
% save('TSOL.mat','TSOL')
%keyboard

%update body orientation
q1=AllData{scaleidx}(4:7,tsol);
q2=AllData{scaleidx}(4:7,tsol+1);
DeltaQ = quatprod(q2,[-q1(1:3);q1(4)]);
qB = quatprod(DeltaQ,Pfull(4:7,end));


%left wing orientation
q1=AllData{scaleidx}(8:11,tsol);
q2=AllData{scaleidx}(8:11,tsol+1);
DeltaQ = quatprod(q2,[-q1(1:3);q1(4)]);
qL = quatprod(DeltaQ,Pfull(8:11,end));


%right wing orientation
q1=AllData{scaleidx}(12:15,tsol);
q2=AllData{scaleidx}(12:15,tsol+1);
DeltaQ = quatprod(q2,[-q1(1:3);q1(4)]);
qR = quatprod(DeltaQ,Pfull(12:15,end));

%transformation
DT = diff(AllData{scaleidx}(1:3,tsol:tsol+1),1,2);

Pnew = [Pfull(1:3,end) + DT
    qB
    qL
    qR];

%keyboard
scaledtime = tt(1):1/scale(scaleidx):tt(2);
tnew = [scaledtime(minidx(scaleidx)) scale(scaleidx)];
    




