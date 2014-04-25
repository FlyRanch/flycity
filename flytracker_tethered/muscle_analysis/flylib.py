import numpy as np
import sys
import json
import os
import random
import scipy
import quantities as pq
import scipy.io
import warnings

stroke_amp_R = 'phi_R'
stroke_amp_L = 'phi_L'
stroke_dev_R = 'theta_R'
stroke_dev_L = 'theta_L'
wing_rot_R = 'eta_R'
wing_rot_L = 'eta_L'

param_file = open('params.json','rb')
params = json.load(param_file)
param_file.close()
#rootpath = params['platform_paths'][sys.platform] + params['root_dir']

class Squadron(object):
    """Controller object to facilitate the groupwise analysis of the fly data"""
    def __init__(self,fly_db):
        self.fly_db = fly_db
        self.flies = [Fly(fly_db[flyn]) for flyn in fly_db.keys()]
    
    def load_kine(self,experiment_name):
        "just load the kine data for each fly"
        for fly in self.flies:
            fly.load_processed_wbkin(experiment_name)
        
class Fly(object):
    """Controler object for fly data, Fly is initialized with a 'fly_record dictionary
    and the object is used to facilitate adding and removing data from this dictionary
    but does does not itself own the dictionary: if Fly is deleted, the dictionary
    can still exist"""
    def __init__(self,fly_db,flynum):
        self.flynum = flynum
        self.fly_record = fly_db[flynum]
        self.param_file = open('params.json','rb')
        self.params = json.load(self.param_file)
        self.param_file.close()
        self.rootpath = self.params['platform_paths'][sys.platform] + self.params['root_dir']
        self.fly_path = self.rootpath + ('Fly%04d/')%(fly_record['flynum'])
    
class Experiment(object):
    """Controller class for an individual experiments init with the fly_record and
    experiment name holds a reference to the experiment record in the fly_record 
    to facilitate operations on those data"""
    
    def __init__(self,fly_record,experiment_name,fly_path):
        self.experiment_name = experiment_name
        self.fly_record = fly_record
        self.exp_record = fly_record['experiments'][experiment_name]
        self.fly_path = fly_path
        
    def load_axon_data(self,filenum = 0):
        """load the axon data from an experiment into the fly_record"""
        axon_file = self.fly_path + self.exp_record['axon_file_names'][filenum]
        self.exp_record['axon_data'] = get_axon_signals(axon_file)
        
    def load_processed_wbkin(self):
        """load the processed (resampled ect) kine data from a cPickle file saved in the 
        primary directory of the fly - gets around the need to load all the axon data
        into memory"""
        import cPickle
        f = self.fly_path + 'lr_blob_expansion_processed_kine.cpkl'
        f = open(f,'rb')
        seqs = cPickle.load(f)
        self.exp_record['wbkin_sequences'] = seqs
        
    def save_processed_wbkin(self):
        """cPickle the processed kine data to be loaded later - helps with memory 
        management""" 
        import cPickle
        f = open(self.fly_path + self.experiment_name + '_processed_wbkin.cpkl','wb')
        cPickle.dump(self.exp_record['wbkin_sequences'],f)
        f.close()   
        
    def load_wbkin_sequences(self):
        """load all the matlab generated kine sequences for an exp"""
        snums = self.exp_record['photron_seq_nums']
        frmtstr = self.exp_record['solution_format_string']
        frmtstr += self.exp_record['kine_filename']
        kinefiles = [self.fly_path + frmtstr%(snum) for snum in snums]
        kines = [load_wbkin_file(kine_filename) for kine_filename in kinefiles]
        #add the file string since if the list of sequences for a fly is re-arranged
        #this key will allow them to be re-sorted in order
        [k.update({'seq_num':snum}) for k,snum in zip(kines,snums)]
        self.exp_record['wbkin_sequences'] = kines
        
    def load_solution_sequences(self):
        """function loads the photron sequences from a given experiment without using
        matlab file"""
        snums = self.exp_record['photron_seq_nums']
        frmtstr = self.exp_record['solution_format_string']
        seqpaths = [self.fly_path + frmtstr%(snum) for snum in snums]
        seqs = [load_solution_files(seqpath) for seqpath in seqpaths]
        self.exp_record['solution_sequences'] = seqs
        
    def get_cam_epochs(self):
        """load and sync the axon and photron data for the capture epochs of a
        given experiment. Loads the axon and photron data if they are not already loaded"""
        fps = pq.Quantity(np.float64(exp['photron_frame_rate_Hz']),'Hz')
        if 'axon_data' not in self.exp_record.keys():
            self.load_axon_data(experiment_name)
        if 'wbkin_sequences' not in self.exp_record.keys(): 
            self.load_wbkin_sequences(experiment_name)
        kine_data = self.exp_record['wbkin_sequences']
        numframes = self.exp_record['wbkin_sequences'][0]['last_track']
        capture_epoch = numframes/fps
        ax_dt = pq.Quantity(self.exp_record['axon_data']['sampling_period'],'s')
        capture_samples = np.ceil(capture_epoch/ax_dt)
        trig_idx = idx_by_thresh(self.exp_record['axon_data']['CamTrig'])
        start_idxs = [x[0] for x in trig_idx]
        cam_epochs = [np.arange(np.int(x),np.int(x)+capture_samples,dtype = np.int) for x in start_idxs]
        def fallback_frame_idx(cam_epoch):
            idx = np.array(np.ceil(np.linspace(cam_epoch[0],cam_epoch[-1],numframes)),dtype = int)
            return idx
        frame_idx_list = list()
        for i,epoch in enumerate(cam_epochs):
            try:
                frame_idxs = self.get_frame_idxs(epoch,self.exp_record['axon_data'])
            except IndexError:
                warnings.warn("problem extracting idxs from camera_sync_signal for"+ \
                "epoch %s using even spaced idx's over the camera epoch instead"%(i))
                frame_idxs = fallback_frame_idx(epoch)
            if not(np.shape(frame_idxs)[0] == numframes):
                import warnings
                warnings.warn("problem extracting idxs from camera_sync_signal for"+ \
                "epoch %s using even spaced idx's over the camera epoch instead"%(i))
                frame_idxs = fallback_frame_idx(epoch)
            frame_idx_list.append(frame_idxs)
        times = exp['axon_data']['times']
        [d.update({'axon_epoch':epoch}) for d,epoch in zip(self.exp_record['wbkin_sequences'],cam_epochs)]
        [d.update({'expan_pol':self.lookup_trial_from_ypos(epoch)}) for d,epoch in zip(self.exp_record['wbkin_sequences'],cam_epochs)]
        [d.update({'axon_idxs':idxs}) for d,idxs in zip(self.exp_record['wbkin_sequences'],frame_idx_list)]
        [d.update({'axon_times':times[idxs]}) for d,idxs in zip(self.exp_record['wbkin_sequences'],frame_idx_list)]
        

class Sequence(object):
    def __init__(self,exp_record,seq_num,fly_path):
        self.exp_record = exp_record
        self.fly_path = fly_path
        self.seq_num = seq_num
        self.seq_record = exp_record[seq_num]
        frmtstr = self.exp_record['solution_format_string']
        self.seq_path = self.fly_path + frmtstr%(self.seq_num)
    
    def load_seqence(self):
        self.seq_record['solution_data'] = load_solution_files(self.seq_path)
        frmtstr = self.exp_record['solution_format_string']
        frmtstr += self.exp_record['kine_filename']
        kine_filename = self.fly_path + frmtstr%(self.seq_num)
        self.seq_record['wbkin_data'] = load_wbkin_file(kine_filename)
        
    def lookup_trial_from_ypos(self,epoch):
        """map the Y position signal to the trial type - given some epoch to average
        over. A future version will be able to figure out what that interval should be
        - but this might be hard to do without loosing generality"""
        if 'axon_data' not in self.exp_record.keys(): self.load_axon_data(experiment_name)
        epoch_ypos = np.mean(self.exp_record['axon_data']['Ypos'][epoch])
        trial_idx = np.argmin(abs(self.exp_record['Ypos_trial_volts']-epoch_ypos))
        trial_val = self.exp_record['Ypos_trial_vals'][trial_idx]
        return trial_val
    
    def calc_seqs_strokeplanes(self):
        """calculate the strokeplane for all the seqences of an experiment
        should move this to the matlab program that generates kine data"""
        if 'solution_sequences' not in self.exp_record.keys(): 
            self.load_solution_sequences()
        seqs = self.exp_record['solution_sequences']
        self.exp_record['strokeplanes'] = [calc_seq_strokeplane(s) for s in seqs]
        
    def resample_strokes(self,seq_num,num_samples = 500):
        """resample the wb into an evenly sampled phase-domain matrix for each
        sequence"""
        wbkin_phases = self.get_kine_phases('lr_blob_expansion',seq_num)
        expmnt = self.fly_record['experiments']['lr_blob_expansion']
        wbkin_keys = [stroke_amp_L,stroke_amp_R,
                     stroke_dev_L,stroke_dev_R,
                     wing_rot_L,wing_rot_R,]
        #resample at these phases
        xi = np.linspace(0,4*np.pi,num_samples)
        #list of sampled wbkin phases
        wbkin_phs_list = wbkin_phases['stroke_phases_kin']
        #list of sampled ephys phases
        ephys_phs_list = wbkin_phases['stroke_phases_axon']
        #list of wbkin times for alignment within the sequence
        times_list = wbkin_phases['stroke_times']
        #list of the kine idxs
        wing_idx_list = wbkin_phases['stroke_kin_idx']
        #list of the ephys idx brackets
        ephys_idx_bracket_list = wbkin_phases['stroke_phys_idx']
        #construct as lists to start
        resampled_strokes = {stroke_amp_L:list(),stroke_amp_R:list(),
                          stroke_dev_L:list(),stroke_dev_R:list(),
                          wing_rot_L:list(),wing_rot_R:list(),
                          'axon_stroke_mtrx':list(),'stroke_idx_in_seq':list(),
                          'stroke_times_in_exp':list()}
        #resample the wb_kine
        from scipy.interpolate import griddata
        for i,idxs in enumerate(wing_idx_list):
            if np.shape(idxs)[0] > 20:
                for kine_key in wbkin_keys:
                    signal = expmnt['wbkin_sequences'][seq_num][kine_key]
                    kine_phase = wbkin_phs_list[i]
                    resampled_signal = griddata(np.unwrap(kine_phase),
                                                signal[idxs],
                                                xi,method = 'cubic')
                    resampled_strokes[kine_key].append(resampled_signal)
                resampled_strokes['stroke_idx_in_seq'].append(i)
                resampled_strokes['stroke_times_in_exp'].append(times_list[i])
        #now add the physiology
        for i in resampled_strokes['stroke_idx_in_seq']:
            x0,x1 = ephys_idx_bracket_list[i]
            ephys_sig = expmnt['axon_data']['AMsysCh1'][x0:x1]
            ephys_phase = ephys_phs_list[i]
            resamp = griddata(ephys_phase,ephys_sig,xi,method = 'cubic')
            resampled_strokes['axon_stroke_mtrx'].append(resamp)
        #convert the lists to arrays
        for key in resampled_strokes.keys():
            resampled_strokes[key] = np.squeeze(np.array(resampled_strokes[key]))
        expmnt['wbkin_sequences'][seq_num]['resampled_strokes'] = resampled_strokes
        
    def get_kine_phases(self,
                        experiment_name,
                        seq_num,
                        mode = 'hilbert',
                        fband = (150,250)):
        """return the time series of the wingstroke phase extracted from the wingbeat
        kine. Send the function the experiment name and the sequence number of intrest
        it should return the data needed to chop up the kine and physiology into 
        wingstrokes and put the data into the phase domain. 'mode' and 'fband' currently
        not used"""
        exp = self.fly_record['experiments'][experiment_name]
        #if 'kine_sequences' not in exp.keys(): self.load_kine_sequences(experiment_name)
        seq = exp['wbkin_sequences'][seq_num]
        from scipy.signal import hilbert
        kine_times = seq['axon_times']
        amp_signal =(seq[stroke_amp_R] + seq[stroke_amp_L])/2
        nan_idx = np.argwhere(~np.isnan(amp_signal))[0][0]
        amp_signal = amp_signal[~np.isnan(amp_signal)]
        kine_times = kine_times[~np.isnan(amp_signal)]
        filt_sig = butter_bandpass_filter(amp_signal,150.,250.,kine_times[1]-kine_times[0],order = 3)
        peaks = scipy.signal.find_peaks_cwt(filt_sig*-1,np.arange(1,20))
        A = np.angle(hilbert(filt_sig))
        A = np.mod(np.unwrap(A),2*np.pi)
        #find the phase of the ventral stroke reversal and re-wrap
        peak_phase = np.mean(A[peaks])
        A2 = np.mod(np.unwrap(A)+peak_phase,2*np.pi)
        idx = np.where(np.diff(A2)<0)[0]+1
        stai = 0
        stpi = len(idx)-2
        #store the data in some lists
        stroke_times = list();stroke_phases = list();stroke_phys_idx = list();
        stroke_kin_idx = list();axon_phases = list()
        #we need the axon data to load the phys data
        if 'axon_data' not in exp.keys(): self.load_axon_data(experiment_name)
        axon_times = exp['axon_data']['times'][seq['axon_epoch']]
        for i1,i2 in zip(idx[stai:stpi],idx[stai+2:stpi+2])[:-3]:
            stroke_times.append(kine_times[i1+nan_idx:i2+nan_idx])
            stroke_phases.append(A2[i1:i2])
            stroke_kin_idx.append(np.arange(i1+nan_idx,i2+nan_idx))
            axon_i1 = np.argwhere(axon_times>=stroke_times[-1][0])[0]+seq['axon_epoch'][0]
            axon_i2 = np.argwhere(axon_times>stroke_times[-1][-1])[0]+seq['axon_epoch'][0]
            stroke_phys_idx.append([axon_i1,axon_i2])
            axon_phases.append(np.linspace(0, 4*np.pi, axon_i2-axon_i1))
        return {'stroke_times':stroke_times,
                'stroke_phases_kin':stroke_phases,
                'stroke_kin_idx':stroke_kin_idx,
                'stroke_phys_idx':stroke_phys_idx,
                'stroke_phases_axon':axon_phases,
                'phot_seq_phases':A2}

def calc_seq_strokeplane(self,seq):
    """calculate the strokeplane from the quaternions of a sequence"""
    q_left = np.squeeze(seq[:,[9,10,11,8]])
    q_right = np.squeeze(seq[:,[13,14,15,12]])
    import transformations as trans
    rot_mats_left = np.array([trans.quaternion_matrix(x.T)[:3,:3] for x in q_left])
    rot_mats_right = np.array([trans.quaternion_matrix(x.T)[:3,:3] for x in q_right])
    l_wingtip = np.dot(rot_mats_left,[0,1,0]).T
    r_wingtip = np.dot(rot_mats_right,[0,-1,0]).T
    from scipy.stats import linregress
    idx = ~np.isnan(l_wingtip[0])
    l_slope = linregress(l_wingtip[0][idx],l_wingtip[2][idx])[0]
    r_slope = linregress(r_wingtip[0][idx],r_wingtip[2][idx])[0]
    return (np.rad2deg(np.arctan(l_slope)),np.rad2deg(np.arctan(r_slope)))
                        
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
    times = blocks[0].segments[0].analogsignals[0].times
    sampling_period = blocks[0].segments[0].analogsignals[0].sampling_period
    times.units = 's'
    sampling_period.units = 's'
    times = np.array(times)
    sampling_period = np.float64(sampling_period)
    signals = [np.array(x) for x in blocks[0].segments[0].analogsignals]
    signals = dict(zip(channel_names,signals))
    signals['times'] = times
    signals['sampling_period'] = sampling_period
    return signals

def load_wbkin_file(kine_filename):
    """load the matlab generated wb kine and append with some info for convenience"""
    kine_data = scipy.io.loadmat(kine_filename)
    ##now to make life easier add the frame numbers to the dictionary
    first_track = np.argwhere(~np.isnan(kine_data['eta_L']))[0][0]+1
    last_track = np.shape(kine_data['eta_R'])[0]
    #sanity check
    def check_first_frame():
        kfn = kine_filename.split('/')[-1]
        flytracks_path = kine_filename.replace(kfn,'flytracks/')
        #get the list of frames
        tracklist = filter(lambda x: x.startswith('fly'), os.listdir(flytracks_path))
        framenums = [np.int(s.strip('fly.mat')) for s in tracklist]
        return np.min(framenums),np.max(framenums)
    assert first_track == check_first_frame()[0],(first_track,check_first_frame()[0])
    assert last_track == check_first_frame()[1],(last_track,check_first_frame()[1])
    kine_data['first_track'] = first_track
    kine_data['last_track'] = last_track
    kine_data['frame_nums'] = np.arange(1,last_track)
    return kine_data

def load_solution_files(sequence_path):
        """convenience function to load the sequence from fly tracks function will cat
        the frame numbers into the first row the kine matrx, also the untracked frames
        before the start of the sequence are padded with NaNs - tries to emulate the same 
        format as the WBkin file"""
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
        
        
def idx_by_thresh(signal,thresh = 0.1):
    idxs = np.squeeze(np.argwhere(signal > thresh))
    split_idxs = np.squeeze(np.argwhere(np.diff(idxs) > 1))
    idx_list = np.split(idxs,split_idxs)
    idx_list = [x[1:] for x in idx_list]
    return idx_list

def butter_bandpass(lowcut, highcut, sampling_period, order=5):
    sampling_frequency = 1.0/sampling_period
    nyq = 0.5 * sampling_frequency
    low = lowcut / nyq
    high = highcut / nyq
    b, a = scipy.signal.butter(order, [low, high], btype='band')
    return b, a

def butter_bandpass_filter(data, lowcut, highcut, sampling_period, order=5):
    b, a = butter_bandpass(lowcut, highcut, sampling_period, order=order)
    y = scipy.signal.filtfilt(b, a, data)
    return y
    
def fit_fourier(strk_mtrx,p_init):
    num_strokes = np.shape(strk_mtrx)[0]
    reshaped = np.squeeze(np.reshape(strk_mtrx,(np.size(strk_mtrx),1)))
    phases = np.linspace(0,2*np.pi*num_strokes,np.size(strk_mtrx))
    y_fit = reshaped[~np.isnan(reshaped)]
    x_fit = phases[~np.isnan(reshaped)]
    from scipy import optimize
    ##speed things up by pre-computing the sin and cosine values
    order = (len(p_init)-1)/2
    n = np.arange(1,order+1)
    onesmat = np.ones((len(n),len(x_fit)))
    phase_mtrx = ((onesmat*x_fit).T*n).T
    cos_mtrx = np.cos(phase_mtrx)
    sin_mtrx = np.sin(phase_mtrx)
    p1,msg = optimize.leastsq(errfunc, p_init[:], args=(cos_mtrx,sin_mtrx,np.rad2deg(y_fit)))
    return p1
    
def fourier_pcomp(p,cos_mtrx,sin_mtrx):
    cp = np.array(p[1:-1:2])[:,np.newaxis]
    sp = np.array(p[2::2])[:,np.newaxis]
    hmtrx = cos_mtrx*cp + sin_mtrx*sp
    return p[0] + np.sum(hmtrx,axis = 0)
    
def fourier(phase,p):
    order = (len(p)-1)/2
    n = np.arange(1,order+1)
    onesmat = np.ones((len(n),len(phase)))
    phase_mtrx = ((onesmat*phase).T*n).T
    cp = np.array(p[1:-1:2])[:,np.newaxis]
    sp = np.array(p[2::2])[:,np.newaxis]
    hmtrx = np.cos(phase_mtrx)*cp + np.sin(phase_mtrx)*sp
    return p[0] + np.sum(hmtrx,axis = 0)

def errfunc(p,cos_mtrx,sin_mtrx,y):
    return fourier_pcomp(p,cos_mtrx,sin_mtrx)-y

