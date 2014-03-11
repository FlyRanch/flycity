clear
clc

% load('flightpathDB_noskip_endMin.mat')
load('flightpathDB_skip10_endMin.mat')

DBseq_all = DBseq;
DBx_all = DBx;
DBy_all = DBy;
DBz_all = DBz;

clear DBseq DBx DBy DBz
% load('flightpathDB_noskip_mid.mat')
load('flightpathDB_skip10_mid.mat')

DBseq_all(end+1:end+size(DBx,2),:) = DBseq;
DBx_all(:,end+1:end+size(DBx,2)) = DBx;
DBy_all(:,end+1:end+size(DBy,2)) = DBy;
DBz_all(:,end+1:end+size(DBz,2)) = DBz;

DBseq = DBseq_all;
DBx = DBx_all;
DBy = DBy_all;
DBz = DBz_all;

save('flightpathDB_skip10_all.mat','DBseq','DBt','DBx','DBy','DBz')