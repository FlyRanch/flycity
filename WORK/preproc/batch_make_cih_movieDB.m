
clear
clc
% warning off
root = cd;

start_path = '/media';
data_path = '/Photron/SEQS';

seq_info_file = 'seq_infoDB.mat';
if exist(seq_info_file) == 2
    load(seq_info_file)
else
    seq_info.date = [];
    seq_info.seq = [];
    seq_info.fps = [];
    seq_info.shutter_inv = [];
    seq_info.N_frames = [];
    seq_info.trigger_frame = [];
    seq_info.im_width = [];
    seq_info.im_height = [];
    seq_info.bits = [];
end

cd(start_path)
discs = dir;
for hd=3:length(discs)
    if exist([discs(hd).name,data_path]) == 7
        cd([discs(hd).name,data_path])
        
        dirs = dir;
        for d=3:length(dirs)
            if dirs(d).isdir==1
                date_now = str2num(dirs(d).name); % dir date
                cd(dirs(d).name)
                subdirs=dir('C001*'); % only first cam
                for s=1:length(subdirs)
                    if subdirs(s).isdir==1
                        seq_now = str2num(subdirs(s).name(end-3:end)); % dir seq
                        cd(subdirs(s).name)

                        infofile = dir('*.cih');
                        if isempty(infofile)==0 

                            % read from infofile
                            cam_info = importdata(infofile(1).name,':',12);
                            fps_rec = cam_info.data(1);
                            cam_info = importdata(infofile(1).name,':',19);
                            save_step = cam_info.data(1);
                            fps = fps_rec/save_step;

                            cam_info = importdata(infofile(1).name,'/',13);
                            shutter_inv = cam_info.data(1);

                            cam_info = importdata(infofile(1).name,':',16);
                            N_frames = cam_info.data(1);
                            frame_start = cam_info.data(2);
                            trigger_frame = 1-frame_start;


                            cam_info = importdata(infofile(1).name,':',20);
                            im_width = cam_info.data(1);
                            im_height = cam_info.data(2);

                            cam_info = importdata(infofile(1).name,':',25);
                            bits = cam_info.data(1);

                            % store in db
                            seq_info.date(end+1,1) = date_now;
                            seq_info.seq(end+1,1) = seq_now;
                            seq_info.fps(end+1,1) = fps;
                            seq_info.shutter_inv(end+1,1) = shutter_inv;
                            seq_info.N_frames(end+1,1) = N_frames;
                            seq_info.trigger_frame(end+1,1) = trigger_frame;
                            seq_info.im_width(end+1,1) = im_width;
                            seq_info.im_height(end+1,1) = im_height;
                            seq_info.bits(end+1,1) = bits;
                        end

                        cd ..

                    end
                end
                cd ..
            end
        end
    end
    cd(start_path)
end
    cd(root)
    save(seq_info_file,'seq_info')