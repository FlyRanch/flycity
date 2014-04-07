import numpy as np
import sys
import json
import os
import random
import scipy
import numpy as np
import quantities as pq
import scipy.io

param_file = open('params.json','rb')
params = json.load(param_file)
param_file.close()
#rootpath = params['platform_paths'][sys.platform] + params['root_dir']

class Fly(object):
    """Controler object for fly data"""
    def __init__(self,fly_record):
        self.flynum = fly_record['flynum']
        self.fly_record = fly_record
        self.param_file = open('params.json','rb')
        self.params = json.load(self.param_file)
        self.param_file.close()
        self.rootpath = self.params['platform_paths'][sys.platform] + self.params['root_dir']
        self.fly_path = self.rootpath + ('Fly%04d/')%(fly_record['flynum'])
        
    def load_photron_sequences(self,experiment_name):
        """function loads the photron sequences from a given experiment"""
        snums = self.fly_record['experiments'][experiment_name]['photron_seq_nums']
        frmtstr = self.fly_record['experiments'][experiment_name]['solution_format_string']
        seqpaths = [self.fly_path + frmtstr%(snum) for snum in snums]
        seqs = [self.load_solution_seq(seqpath) for seqpath in seqpaths]
        self.fly_record['experiments'][experiment_name]['photron_sequences'] = seqs
        
    def load_solution_seq(self,sequence_path):
        """convenience function to load the sequence from fly tracks function will concatinate
        the frame numbers into the first row the kine matrxi, also the untracked before the
        start of the sequence are padded with NaNs"""
        flytracks_path = sequence_path + 'flytracks/'
        #get the list of frames
        tracklist = filter(lambda x: x.startswith('fly'), os.listdir(flytracks_path))
        framenums = [np.int(s.strip('fly.mat')) for s in tracklist]
        #sort the files by frame number
        tracklist = [x for (y,x) in sorted(zip(framenums,tracklist))]
        tracklist = [flytracks_path + x for x in tracklist]
        #load the data from the .mat files
        mats = [scipy.io.loadmat(trk)['xh'].copy() for trk in tracklist]
        mats = np.array([np.squeeze(np.pad(np.squeeze(mat),(0,15-np.shape(mat)[0]),'constant')) for mat in mats])
        #since the frames are sorted - now sort the frame numbers
        framenums = np.squeeze(sorted(framenums))
        #cat the frame numbers onto the kine data
        seq = np.concatenate((framenums[:,np.newaxis],mats),axis = 1)
        #now pad the matrix with NaNs so that it starts at frame 0
        start_frame = seq[0,0]
        pad_seq = np.pad(seq,((start_frame-1,0),(0,0)),'constant')
        pad_seq[:start_frame-1,0] = np.arange(1,start_frame)
        pad_seq[:start_frame-1,1:] = np.NAN
        return pad_seq

    def load_axon_data(self,experiment_name,filenum = 0):
        """load the axon data from an experiment into the fly_record"""
        axon_file = self.fly_path + self.fly_record['experiments'][experiment_name]['axon_file_names'][filenum]
        self.fly_record['experiments'][experiment_name]['axon_data'] = get_axon_signals(axon_file) 
    
    def get_cam_epochs(self,experiment_name,numframes = 1365,fps = pq.Quantity(6000.0,'Hz')):
        """load and sync the axon and photron data for the capture epochs of a
        given experiment. Loads the axon data if it is not already loaded"""
        if 'axon_data' in self.fly_record['experiments'][experiment_name].keys():
            axondata = self.fly_record['experiments'][experiment_name]['axon_data']
        else:
            self.load_axon_data(experiment_name)
            axondata = self.fly_record['experiments'][experiment_name]['axon_data']
        capture_epoch = numframes/fps
        ax_dt = axondata['CamTrig'].sampling_period
        ax_dt.units = 's'
        capture_samples = np.ceil(capture_epoch/ax_dt)
        trig_idx = idx_by_thresh(axondata['CamTrig'])
        start_idxs = [x[0] for x in trig_idx]
        cam_epochs = [np.arange(np.int(x),np.int(x)+capture_samples,dtype = np.int) for x in start_idxs]
        frame_idxs = [self.get_frame_idxs(epoch,axondata) for epoch in cam_epochs]
        times = axondata['CamTrig'].times
        self.fly_record['experiments'][experiment_name]['cam_epochs'] = cam_epochs
        self.fly_record['experiments'][experiment_name]['frame_idxs'] = frame_idxs
        self.fly_record['experiments'][experiment_name]['frame_times'] = [times[idx] for idx in frame_idxs]
    
    def get_frame_idxs(self,cam_epoch,axondata):
        """exctract the sync pulse from the camera epochs"""
        frame_idxs = [x[0]+cam_epoch[0] for x in idx_by_thresh(axondata['CamSync'][cam_epoch]*-1,-3.5)]
        frame_idxs[0] -=1
        return frame_idxs


def get_axon_signals(filename):
    from neo.io.axonio import AxonIO
    reader = AxonIO(filename=filename)
    blocks = reader.read()
    header = reader.read_header()
    channel_names = [info['ADCChNames'] for info in header['listADCInfo']]
    channel_units = [info['ADCChUnits'] for info in header['listADCInfo']]
    signals = dict(zip(channel_names,blocks[0].segments[0].analogsignals))
    return signals

def idx_by_thresh(signal,thresh = 0.1):
    idxs = np.squeeze(np.argwhere(signal > thresh))
    split_idxs = np.squeeze(np.argwhere(np.diff(idxs) > 1))
    idx_list = np.split(idxs,split_idxs)
    idx_list = [x[1:] for x in idx_list]
    return idx_list


