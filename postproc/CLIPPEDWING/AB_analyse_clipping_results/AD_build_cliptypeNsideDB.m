clear
clc

clip_typeNside.info.type{1,1}='distal_tip';
clip_typeNside.info.type{2,1}='trailing_edge';

clip_typeNside.info.side{1,1}='left';
clip_typeNside.info.side{2,1}='right';

T=readtable('clip_typeNside.csv','ReadVariableNames',false)

for i=1:size(T,1)
    
    date_now = T{i,1};
    date_now = cell2mat(date_now);
    date_now = str2double([date_now(1:4) date_now(6:7) date_now(9:10)]);
    
    clip_typeNside.date(i,1) = date_now;
end

% clip_typeNside.date = table2struct(T(:,1));
clip_typeNside.seq = table2array(T(:,2));
clip_typeNside.clip_type = table2array(T(:,3));
clip_typeNside.clip_side = table2array(T(:,4));
    
save('clip_typeNsideDB.mat','clip_typeNside')