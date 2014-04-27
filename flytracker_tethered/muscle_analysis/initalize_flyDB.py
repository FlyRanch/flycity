# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

import sys
import json
import os
import numpy as np
import h5py

def main():
    import cPickle
    fn = 'fly_db_init.cpkl'
    f = open(fn,'wb')
    cPickle.dump(fly_db,f)
    f.close()
    
    
fly_db = h5py.File("flydb.hdf5", "w")
#############################################################################################fly_record = dict()
flynum = str(111)
fly_db.create_group(flynum)
fly_record =fly_db[flynum]
fly_record['flynum'] = flynum
fly_record.create_group('experiments')
fly_record['experiments'].create_group('lr_blob_expansion')
#fly_record['experiments'].create_group('b1_azm_expansion_tuning')
fly_record['experiments']['lr_blob_expansion']['photron_seq_nums'] = [1,2,3,4,5,6]
fly_record['experiments']['lr_blob_expansion']['axon_file_names'] = ['fly01_lr_blob_expansion_14401000.abf']
fly_record['experiments']['lr_blob_expansion']['photron_date_string'] = ['20140401']
fly_record['experiments']['lr_blob_expansion']['kine_filename'] = ['WBkin.mat']
fly_record['experiments']['lr_blob_expansion']['solution_format_string'] = ['20140401_S%04d/']
fly_record['experiments']['lr_blob_expansion']['photron_frame_rate_Hz'] = 6000
fly_record['experiments']['lr_blob_expansion']['Ypos_trial_volts'] = np.linspace(1,10,12)
fly_record['experiments']['lr_blob_expansion']['Ypos_trial_vals'] = np.concatenate(([np.nan],np.arange(0,12)*30))
fly_record['experiments']['lr_blob_expansion'].create_group('sequences')
#fly_record['experiments']['b1_azm_expansion_tuning']['axon_file_names'] = ['fly01_b1tuning_14401001.abf']
#############################################################################################
flynum = str(112)
fly_db.create_group(flynum)
fly_record =fly_db[flynum]
fly_record['flynum'] = flynum
fly_record.create_group('experiments')
fly_record['experiments'].create_group('lr_blob_expansion')
#fly_record['experiments'].create_group('b1_azm_expansion_tuning')
fly_record['experiments']['lr_blob_expansion']['photron_seq_nums'] = [7,8,9,10,11,12]
fly_record['experiments']['lr_blob_expansion']['axon_file_names'] = ['fly02_lr_blob_expansion_14401002.abf']
fly_record['experiments']['lr_blob_expansion']['photron_date_string'] = ['20140401']
fly_record['experiments']['lr_blob_expansion']['kine_filename'] = ['WBkin.mat']
fly_record['experiments']['lr_blob_expansion']['solution_format_string'] = ['20140401_S%04d/']
fly_record['experiments']['lr_blob_expansion']['photron_frame_rate_Hz'] = 6000
fly_record['experiments']['lr_blob_expansion']['Ypos_trial_volts'] = np.linspace(1,10,12)
fly_record['experiments']['lr_blob_expansion']['Ypos_trial_vals'] = np.concatenate(([np.nan],np.arange(0,12)*30))
fly_record['experiments']['lr_blob_expansion'].create_group('sequences')
#fly_record['experiments']['b1_azm_expansion_tuning']['axon_file_names'] = ['fly02_b1tuning_14401003.abf']
#############################################################################################
flynum = str(114)
fly_db.create_group(flynum)
fly_record =fly_db[flynum]
fly_record['flynum'] = flynum
fly_record.create_group('experiments')
fly_record['experiments'].create_group('lr_blob_expansion')
#fly_record['experiments'].create_group('b1_azm_expansion_tuning')
fly_record['experiments']['lr_blob_expansion']['photron_seq_nums'] = [13,14,15,16,17,18]
fly_record['experiments']['lr_blob_expansion']['axon_file_names'] = ['fly04_lr_blob_expansion_14401012.abf']
fly_record['experiments']['lr_blob_expansion']['photron_date_string'] = ['20140401']
fly_record['experiments']['lr_blob_expansion']['kine_filename'] = ['WBkin.mat']
fly_record['experiments']['lr_blob_expansion']['solution_format_string'] = ['20140401_S%04d/']
fly_record['experiments']['lr_blob_expansion']['photron_frame_rate_Hz'] = 6000
fly_record['experiments']['lr_blob_expansion']['Ypos_trial_volts'] = np.linspace(1,10,12)
fly_record['experiments']['lr_blob_expansion']['Ypos_trial_vals'] = np.concatenate(([np.nan],np.arange(0,12)*30))
fly_record['experiments']['lr_blob_expansion'].create_group('sequences')
#fly_record['experiments']['b1_azm_expansion_tuning']['axon_file_names'] = ['fly01_b1tuning_14401013.abf']
#############################################################################################
flynum = str(115)
fly_db.create_group(flynum)
fly_record =fly_db[flynum]
fly_record['flynum'] = flynum
fly_record.create_group('experiments')
fly_record['experiments'].create_group('lr_blob_expansion')
#fly_record['experiments'].create_group('b1_azm_expansion_tuning')
fly_record['experiments']['lr_blob_expansion']['photron_seq_nums'] = [1,2,3,4,5,6]
fly_record['experiments']['lr_blob_expansion']['axon_file_names'] = ['fly01_lr_blob_expansion_14402001.abf']
fly_record['experiments']['lr_blob_expansion']['photron_date_string'] = ['20140402']
fly_record['experiments']['lr_blob_expansion']['kine_filename'] = ['WBkin.mat']
fly_record['experiments']['lr_blob_expansion']['solution_format_string'] = ['20140402_S%04d/']
fly_record['experiments']['lr_blob_expansion']['photron_frame_rate_Hz'] = 6000
fly_record['experiments']['lr_blob_expansion']['Ypos_trial_volts'] = np.linspace(1,10,12)
fly_record['experiments']['lr_blob_expansion']['Ypos_trial_vals'] = np.concatenate(([np.nan],np.arange(0,12)*30))
fly_record['experiments']['lr_blob_expansion'].create_group('sequences')
#fly_record['experiments']['b1_azm_expansion_tuning']['axon_file_names'] = ['fly01_b1tuning_14402002.abf']

#if __name__ == '__main__':
#s    main()
