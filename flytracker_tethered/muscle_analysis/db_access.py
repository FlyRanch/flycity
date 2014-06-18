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

def get_db():
    fly_db = h5py.File('/Volumes/FlyDataB/FlyDB/flydb.hdf5','a')
    return fly_db

starfield_pattern_names_6_0_2014  = ['equator_000.mat',
                            'equator_030.mat',
                            'equator_060.mat',
                            'equator_090.mat',
                            'equator_120.mat',
                            'equator_150.mat',
                            'equator_180.mat',
                            'equator_210.mat',
                            'equator_240.mat',
                            'equator_270.mat',
                            'equator_300.mat',
                            'equator_330.mat',
                            'coromeridian_030.mat',
                            'coromeridian_060.mat',
                            'coromeridian_090.mat',
                            'coromeridian_120.mat',
                            'coromeridian_150.mat',
                            'coromeridian_210.mat',
                            'coromeridian_240.mat',
                            'coromeridian_270.mat',
                            'coromeridian_300.mat',
                            'coromeridian_330.mat',
                            'sagimeridian_030.mat',
                            'sagimeridian_060.mat',
                            'sagimeridian_120.mat',
                            'sagimeridian_150.mat',
                            'sagimeridian_210.mat',
                            'sagimeridian_240.mat',
                            'sagimeridian_300.mat',
                            'sagimeridian_330.mat']
    
def init_db():
    fly_db = h5py.File("/Vlumes/FlyDataB/FlyDB/flydb.hdf5", "w")
    #############################################################################################fly_record = dict()
    flynum = str(111)
    fly_db.create_group(flynum)
    fly_record =fly_db[flynum]
    fly_record['flynum'] = flynum
    fly_record.create_group('experiments')
    fly_record['experiments'].create_group('lr_blob_expansion')
    fly_record['experiments']['lr_blob_expansion']['photron_seq_nums'] = [1,2,3,4,5,6]
    fly_record['experiments']['lr_blob_expansion']['axon_file_names'] = ['fly01_lr_blob_expansion_14401000.abf']
    fly_record['experiments']['lr_blob_expansion']['photron_date_string'] = ['20140401']
    fly_record['experiments']['lr_blob_expansion']['kine_filename'] = ['WBkin.mat']
    fly_record['experiments']['lr_blob_expansion']['solution_format_string'] = ['20140401_S%04d/']
    fly_record['experiments']['lr_blob_expansion']['photron_frame_rate_Hz'] = 6000
    fly_record['experiments']['lr_blob_expansion']['Ypos_trial_volts'] = np.linspace(1,10,12)
    fly_record['experiments']['lr_blob_expansion']['Ypos_trial_vals'] = np.concatenate(([np.nan],np.arange(0,12)*30))
    fly_record['experiments']['lr_blob_expansion']['AMsysCh1_ID'] = 'b1'
    fly_record['experiments']['lr_blob_expansion']['AMsysCh1_side'] = 'l'
    fly_record['experiments']['lr_blob_expansion'].create_group('sequences')
    #############################################################################################
    flynum = str(112)
    fly_db.create_group(flynum)
    fly_record =fly_db[flynum]
    fly_record['flynum'] = flynum
    fly_record.create_group('experiments')
    fly_record['experiments'].create_group('lr_blob_expansion')
    fly_record['experiments']['lr_blob_expansion']['photron_seq_nums'] = [7,8,9,10,11,12]
    fly_record['experiments']['lr_blob_expansion']['axon_file_names'] = ['fly02_lr_blob_expansion_14401002.abf']
    fly_record['experiments']['lr_blob_expansion']['photron_date_string'] = ['20140401']
    fly_record['experiments']['lr_blob_expansion']['kine_filename'] = ['WBkin.mat']
    fly_record['experiments']['lr_blob_expansion']['solution_format_string'] = ['20140401_S%04d/']
    fly_record['experiments']['lr_blob_expansion']['photron_frame_rate_Hz'] = 6000
    fly_record['experiments']['lr_blob_expansion']['Ypos_trial_volts'] = np.linspace(1,10,12)
    fly_record['experiments']['lr_blob_expansion']['Ypos_trial_vals'] = np.concatenate(([np.nan],np.arange(0,12)*30))
    fly_record['experiments']['lr_blob_expansion']['AMsysCh1_ID'] = 'b1'
    fly_record['experiments']['lr_blob_expansion']['AMsysCh1_side'] = 'l'
    fly_record['experiments']['lr_blob_expansion'].create_group('sequences')
    #############################################################################################
    flynum = str(114)
    fly_db.create_group(flynum)
    fly_record =fly_db[flynum]
    fly_record['flynum'] = flynum
    fly_record.create_group('experiments')
    fly_record['experiments'].create_group('lr_blob_expansion')
    fly_record['experiments']['lr_blob_expansion']['photron_seq_nums'] = [13,14,15,16,17,18]
    fly_record['experiments']['lr_blob_expansion']['axon_file_names'] = ['fly04_lr_blob_expansion_14401012.abf']
    fly_record['experiments']['lr_blob_expansion']['photron_date_string'] = ['20140401']
    fly_record['experiments']['lr_blob_expansion']['kine_filename'] = ['WBkin.mat']
    fly_record['experiments']['lr_blob_expansion']['solution_format_string'] = ['20140401_S%04d/']
    fly_record['experiments']['lr_blob_expansion']['photron_frame_rate_Hz'] = 6000
    fly_record['experiments']['lr_blob_expansion']['Ypos_trial_volts'] = np.linspace(1,10,12)
    fly_record['experiments']['lr_blob_expansion']['Ypos_trial_vals'] = np.concatenate(([np.nan],np.arange(0,12)*30))
    fly_record['experiments']['lr_blob_expansion']['AMsysCh1_ID'] = 'b1'
    fly_record['experiments']['lr_blob_expansion']['AMsysCh1_side'] = 'l'
    fly_record['experiments']['lr_blob_expansion'].create_group('sequences')
    #############################################################################################
    flynum = str(115)
    fly_db.create_group(flynum)
    fly_record =fly_db[flynum]
    fly_record['flynum'] = flynum
    fly_record.create_group('experiments')
    fly_record['experiments'].create_group('lr_blob_expansion')
    fly_record['experiments']['lr_blob_expansion']['photron_seq_nums'] = [1,2,3,4,5,6]
    fly_record['experiments']['lr_blob_expansion']['axon_file_names'] = ['fly01_lr_blob_expansion_14402001.abf']
    fly_record['experiments']['lr_blob_expansion']['photron_date_string'] = ['20140402']
    fly_record['experiments']['lr_blob_expansion']['kine_filename'] = ['WBkin.mat']
    fly_record['experiments']['lr_blob_expansion']['solution_format_string'] = ['20140402_S%04d/']
    fly_record['experiments']['lr_blob_expansion']['photron_frame_rate_Hz'] = 6000
    fly_record['experiments']['lr_blob_expansion']['Ypos_trial_volts'] = np.linspace(1,10,12)
    fly_record['experiments']['lr_blob_expansion']['Ypos_trial_vals'] = np.concatenate(([np.nan],np.arange(0,12)*30))
    fly_record['experiments']['lr_blob_expansion']['AMsysCh1_ID'] = 'b1'
    fly_record['experiments']['lr_blob_expansion']['AMsysCh1_side'] = 'l'
    fly_record['experiments']['lr_blob_expansion'].create_group('sequences')
    #############################################################################################
    flynum = str(116)
    fly_db.create_group(flynum)
    fly_record =fly_db[flynum]
    fly_record['flynum'] = flynum
    fly_record.create_group('experiments')
    fly_record['experiments'].create_group('lr_blob_expansion')
    fly_record['experiments']['lr_blob_expansion']['photron_seq_nums'] = [1,2,3,4,5,6]
    fly_record['experiments']['lr_blob_expansion']['axon_file_names'] = ['fly01_lr_blob_expansion_14410000.abf']
    fly_record['experiments']['lr_blob_expansion']['photron_date_string'] = ['20140410']
    fly_record['experiments']['lr_blob_expansion']['kine_filename'] = ['WBkin.mat']
    fly_record['experiments']['lr_blob_expansion']['solution_format_string'] = ['20140410_S%04d/']
    fly_record['experiments']['lr_blob_expansion']['photron_frame_rate_Hz'] = 6000
    fly_record['experiments']['lr_blob_expansion']['Ypos_trial_volts'] = np.linspace(1,10,12)
    fly_record['experiments']['lr_blob_expansion']['Ypos_trial_vals'] = np.concatenate(([np.nan],np.arange(0,12)*30))
    fly_record['experiments']['lr_blob_expansion']['AMsysCh1_ID'] = 'b1'
    fly_record['experiments']['lr_blob_expansion']['AMsysCh1_side'] = 'l'
    fly_record['experiments']['lr_blob_expansion'].create_group('sequences')
    #############################################################################################
    flynum = str(117)
    fly_db.create_group(flynum)
    fly_record =fly_db[flynum]
    fly_record['flynum'] = flynum
    fly_record.create_group('experiments')
    fly_record['experiments'].create_group('lr_blob_expansion')
    #fly_record['experiments'].create_group('b1_azm_expansion_tuning')
    fly_record['experiments']['lr_blob_expansion']['photron_seq_nums'] = [7,8,9,10,11,12]
    fly_record['experiments']['lr_blob_expansion']['axon_file_names'] = ['fly02_lr_blob_expansion_14410002.abf']
    fly_record['experiments']['lr_blob_expansion']['photron_date_string'] = ['20140410']
    fly_record['experiments']['lr_blob_expansion']['kine_filename'] = ['WBkin.mat']
    fly_record['experiments']['lr_blob_expansion']['solution_format_string'] = ['20140410_S%04d/']
    fly_record['experiments']['lr_blob_expansion']['photron_frame_rate_Hz'] = 6000
    fly_record['experiments']['lr_blob_expansion']['Ypos_trial_volts'] = np.linspace(1,10,12)
    fly_record['experiments']['lr_blob_expansion']['Ypos_trial_vals'] = np.concatenate(([np.nan],np.arange(0,12)*30))
    fly_record['experiments']['lr_blob_expansion']['AMsysCh1_ID'] = 'b1'
    fly_record['experiments']['lr_blob_expansion']['AMsysCh1_side'] = 'l'
    fly_record['experiments']['lr_blob_expansion'].create_group('sequences')
    #############################################################################################
    flynum = str(118)
    fly_db.create_group(flynum)
    fly_record =fly_db[flynum]
    fly_record['flynum'] = flynum
    fly_record.create_group('experiments')
    fly_record['experiments'].create_group('lr_blob_expansion')
    #fly_record['experiments'].create_group('b1_azm_expansion_tuning')
    fly_record['experiments']['lr_blob_expansion']['photron_seq_nums'] = [13,14,15,16,17]
    fly_record['experiments']['lr_blob_expansion']['axon_file_names'] = ['fly03_lr_blob_expansion_14410003.abf']
    fly_record['experiments']['lr_blob_expansion']['photron_date_string'] = ['20140410']
    fly_record['experiments']['lr_blob_expansion']['kine_filename'] = ['WBkin.mat']
    fly_record['experiments']['lr_blob_expansion']['solution_format_string'] = ['20140410_S%04d/']
    fly_record['experiments']['lr_blob_expansion']['photron_frame_rate_Hz'] = 6000
    fly_record['experiments']['lr_blob_expansion']['Ypos_trial_volts'] = np.linspace(1,10,12)
    fly_record['experiments']['lr_blob_expansion']['Ypos_trial_vals'] = np.concatenate(([np.nan],np.arange(0,12)*30))
    fly_record['experiments']['lr_blob_expansion']['AMsysCh1_ID'] = 'b1'
    fly_record['experiments']['lr_blob_expansion']['AMsysCh1_side'] = 'l'
    fly_record['experiments']['lr_blob_expansion'].create_group('sequences')
    #############################################################################################
    flynum = str(122)
    fly_db.create_group(flynum)
    fly_record =fly_db[flynum]
    fly_record['flynum'] = flynum
    fly_record.create_group('experiments')
    fly_record['experiments'].create_group('lr_blob_expansion')
    #fly_record['experiments'].create_group('b1_azm_expansion_tuning')
    fly_record['experiments']['lr_blob_expansion']['photron_seq_nums'] = [1,2,3,4,5,6]
    fly_record['experiments']['lr_blob_expansion']['axon_file_names'] = ['fly01_lr_blob_expansion_14428000.abf']
    fly_record['experiments']['lr_blob_expansion']['photron_date_string'] = ['20140428']
    fly_record['experiments']['lr_blob_expansion']['kine_filename'] = ['WBkin.mat']
    fly_record['experiments']['lr_blob_expansion']['solution_format_string'] = ['20140428_S%04d/']
    fly_record['experiments']['lr_blob_expansion']['photron_frame_rate_Hz'] = 6000
    fly_record['experiments']['lr_blob_expansion']['Ypos_trial_volts'] = np.linspace(1,10,12)
    fly_record['experiments']['lr_blob_expansion']['Ypos_trial_vals'] = np.concatenate(([np.nan],np.arange(0,12)*30))
    fly_record['experiments']['lr_blob_expansion']['AMsysCh1_ID'] = 'b1'
    fly_record['experiments']['lr_blob_expansion']['AMsysCh1_side'] = 'r'
    fly_record['experiments']['lr_blob_expansion'].create_group('sequences')
    #############################################################################################
    flynum = str(123)
    fly_db.create_group(flynum)
    fly_record =fly_db[flynum]
    fly_record['flynum'] = flynum
    fly_record.create_group('experiments')
    fly_record['experiments'].create_group('lr_blob_expansion')
    #fly_record['experiments'].create_group('b1_azm_expansion_tuning')
    fly_record['experiments']['lr_blob_expansion']['photron_seq_nums'] = [7,8,9,10,11]
    fly_record['experiments']['lr_blob_expansion']['axon_file_names'] = ['fly02_lr_blob_expansion_14428001.abf']
    fly_record['experiments']['lr_blob_expansion']['photron_date_string'] = ['20140428']
    fly_record['experiments']['lr_blob_expansion']['kine_filename'] = ['WBkin.mat']
    fly_record['experiments']['lr_blob_expansion']['solution_format_string'] = ['20140428_S%04d/']
    fly_record['experiments']['lr_blob_expansion']['photron_frame_rate_Hz'] = 6000
    fly_record['experiments']['lr_blob_expansion']['Ypos_trial_volts'] = np.linspace(1,10,12)
    fly_record['experiments']['lr_blob_expansion']['Ypos_trial_vals'] = np.concatenate(([np.nan],np.arange(0,12)*30))
    fly_record['experiments']['lr_blob_expansion']['AMsysCh1_ID'] = 'b1'
    fly_record['experiments']['lr_blob_expansion']['AMsysCh1_side'] = 'r'
    fly_record['experiments']['lr_blob_expansion'].create_group('sequences')
    #############################################################################################
    flynum = str(124)
    fly_db.create_group(flynum)
    fly_record =fly_db[flynum]
    fly_record['flynum'] = flynum
    fly_record.create_group('experiments')
    fly_record['experiments'].create_group('lr_blob_expansion')
    #fly_record['experiments'].create_group('b1_azm_expansion_tuning')
    fly_record['experiments']['lr_blob_expansion']['photron_seq_nums'] = [2,3,4,5,6]
    fly_record['experiments']['lr_blob_expansion']['axon_file_names'] = ['fly01_lr_blob_expansion_14429001.abf']
    fly_record['experiments']['lr_blob_expansion']['photron_date_string'] = ['20140429']
    fly_record['experiments']['lr_blob_expansion']['kine_filename'] = ['WBkin.mat']
    fly_record['experiments']['lr_blob_expansion']['solution_format_string'] = ['20140429_S%04d/']
    fly_record['experiments']['lr_blob_expansion']['photron_frame_rate_Hz'] = 6000
    fly_record['experiments']['lr_blob_expansion']['Ypos_trial_volts'] = np.linspace(1,10,12)
    fly_record['experiments']['lr_blob_expansion']['Ypos_trial_vals'] = np.concatenate(([np.nan],np.arange(0,12)*30))
    fly_record['experiments']['lr_blob_expansion']['AMsysCh1_ID'] = 'b1'
    fly_record['experiments']['lr_blob_expansion']['AMsysCh1_side'] = 'r'
    fly_record['experiments']['lr_blob_expansion'].create_group('sequences')
    #############################################################################################
    flynum = str(125)
    fly_db.create_group(flynum)
    fly_record =fly_db[flynum]
    fly_record['flynum'] = flynum
    fly_record.create_group('experiments')
    fly_record['experiments'].create_group('lr_blob_expansion')
    #fly_record['experiments'].create_group('b1_azm_expansion_tuning')
    fly_record['experiments']['lr_blob_expansion']['photron_seq_nums'] = [7,8,9,10,11,12]
    fly_record['experiments']['lr_blob_expansion']['axon_file_names'] = ['fly02_lr_blob_expansion_14429004.abf']
    fly_record['experiments']['lr_blob_expansion']['photron_date_string'] = ['20140429']
    fly_record['experiments']['lr_blob_expansion']['kine_filename'] = ['WBkin.mat']
    fly_record['experiments']['lr_blob_expansion']['solution_format_string'] = ['20140429_S%04d/']
    fly_record['experiments']['lr_blob_expansion']['photron_frame_rate_Hz'] = 6000
    fly_record['experiments']['lr_blob_expansion']['Ypos_trial_volts'] = np.linspace(1,10,12)
    fly_record['experiments']['lr_blob_expansion']['Ypos_trial_vals'] = np.concatenate(([np.nan],np.arange(0,12)*30))
    fly_record['experiments']['lr_blob_expansion']['AMsysCh1_ID'] = 'b1'
    fly_record['experiments']['lr_blob_expansion']['AMsysCh1_side'] = 'r'
    fly_record['experiments']['lr_blob_expansion'].create_group('sequences')
    #############################################################################################
    #############################################################################################
    #############################################################################################
    #############################################################################################

    flynum = str(130)
    fly_db.create_group(flynum)
    fly_record =fly_db[flynum]
    fly_record['flynum'] = flynum
    fly_record.create_group('experiments')
    fly_record['experiments'].create_group('lr_blob_expansion')
    #fly_record['experiments'].create_group('b1_azm_expansion_tuning')
    fly_record['experiments']['lr_blob_expansion']['photron_seq_nums'] = [1,2,3,4,5,6]
    fly_record['experiments']['lr_blob_expansion']['axon_file_names'] = ['fly02_lr_blob_expansion_14429004.abf']
    fly_record['experiments']['lr_blob_expansion']['photron_date_string'] = ['20140429']
    fly_record['experiments']['lr_blob_expansion']['kine_filename'] = ['WBkin.mat']
    fly_record['experiments']['lr_blob_expansion']['solution_format_string'] = ['20140429_S%04d/']
    fly_record['experiments']['lr_blob_expansion']['photron_frame_rate_Hz'] = 6000
    fly_record['experiments']['lr_blob_expansion']['Ypos_trial_volts'] = np.linspace(1,10,12)
    fly_record['experiments']['lr_blob_expansion']['Ypos_trial_vals'] = np.concatenate(([np.nan],np.arange(0,12)*30))
    fly_record['experiments']['lr_blob_expansion']['AMsysCh1_ID'] = 'i1'
    fly_record['experiments']['lr_blob_expansion']['AMsysCh1_side'] = 'r'
    fly_record['experiments']['lr_blob_expansion'].create_group('sequences')

    #############################################################################################
    #############################################################################################


    flynum = str(151)
    fly_db.create_group(flynum)
    fly_record =fly_db[flynum]
    fly_record['flynum'] = flynum
    fly_record.create_group('experiments')
    fly_record['experiments'].create_group('img_starfield_t2_rep1')
    #fly_record['experiments'].create_group('b1_azm_expansion_tuning')
    fly_record['experiments']['img_starfield_t2_rep1']['axon_file_names'] = ['fly0151_rotating_starfield_imaging_T2_trial_1_14529002.abf']
    fly_record['experiments']['img_starfield_t2_rep1']['tiff_file_names'] = ['/trial1/trial1_MMStack.ome.tif']
    fly_record['experiments']['img_starfield_t2_rep1']['sequence_pattern_names'] = starfield_pattern_names_6_0_2014
    fly_record['experiments']['img_starfield_t2_rep1'].create_group('sequences')

    flynum = str(153)
    fly_db.create_group(flynum)
    fly_record =fly_db[flynum]
    fly_record['flynum'] = flynum
    fly_record.create_group('experiments')
    fly_record['experiments'].create_group('img_starfield_t2_rep1')
    #fly_record['experiments'].create_group('b1_azm_expansion_tuning')
    fly_record['experiments']['img_starfield_t2_rep1']['axon_file_names'] = ['fly0153_rotating_starfield_imaging_T2_trial_1_14530005.abf']
    fly_record['experiments']['img_starfield_t2_rep1']['tiff_file_names'] = ['/T2_trial1/T2_trial1_MMStack.ome.tif']
    fly_record['experiments']['img_starfield_t2_rep1']['sequence_pattern_names'] = starfield_pattern_names_6_0_2014
    fly_record['experiments']['img_starfield_t2_rep1'].create_group('sequences')

    flynum = str(154)
    fly_db.create_group(flynum)
    fly_record =fly_db[flynum]
    fly_record['flynum'] = flynum
    fly_record.create_group('experiments')
    fly_record['experiments'].create_group('img_starfield_t2_rep1')
    #fly_record['experiments'].create_group('b1_azm_expansion_tuning')
    fly_record['experiments']['img_starfield_t2_rep1']['axon_file_names'] = ['fly0154_rotating_starfield_imaging_T2_trial_1_14530007.abf']
    fly_record['experiments']['img_starfield_t2_rep1']['tiff_file_names'] = ['/T2_trial1/T2_trial1_MMStack.ome.tif']
    fly_record['experiments']['img_starfield_t2_rep1']['sequence_pattern_names'] = starfield_pattern_names_6_0_2014
    fly_record['experiments']['img_starfield_t2_rep1'].create_group('sequences')


    flynum = str(155)
    fly_db.create_group(flynum)
    fly_record =fly_db[flynum]
    fly_record['flynum'] = flynum
    fly_record.create_group('experiments')
    fly_record['experiments'].create_group('img_starfield_t2_rep1')
    #fly_record['experiments'].create_group('b1_azm_expansion_tuning')
    fly_record['experiments']['img_starfield_t2_rep1']['axon_file_names'] = ['fly0155_rotating_starfield_imaging_T2_trial_1_14530009.abf']
    fly_record['experiments']['img_starfield_t2_rep1']['tiff_file_names'] = ['/T2_trial1/T2_trial1_MMStack.ome.tif']
    fly_record['experiments']['img_starfield_t2_rep1']['sequence_pattern_names'] = starfield_pattern_names_6_0_2014
    fly_record['experiments']['img_starfield_t2_rep1'].create_group('sequences')


    flynum = str(156)
    fly_db.create_group(flynum)
    fly_record =fly_db[flynum]
    fly_record['flynum'] = flynum
    fly_record.create_group('experiments')
    fly_record['experiments'].create_group('img_starfield_t2_rep1')
    #fly_record['experiments'].create_group('b1_azm_expansion_tuning')
    fly_record['experiments']['img_starfield_t2_rep1']['axon_file_names'] = ['fly0156_rotating_starfield_imaging_T2_trial_1_14530011.abf']
    fly_record['experiments']['img_starfield_t2_rep1']['tiff_file_names'] = ['/T2_trial1/T2_trial1_MMStack.ome.tif']
    fly_record['experiments']['img_starfield_t2_rep1']['sequence_pattern_names'] = starfield_pattern_names_6_0_2014
    fly_record['experiments']['img_starfield_t2_rep1'].create_group('sequences')

    flynum = str(157)
    fly_db.create_group(flynum)
    fly_record =fly_db[flynum]
    fly_record['flynum'] = flynum
    fly_record.create_group('experiments')
    fly_record['experiments'].create_group('img_starfield_t2_rep1')
    #fly_record['experiments'].create_group('b1_azm_expansion_tuning')
    fly_record['experiments']['img_starfield_t2_rep1']['axon_file_names'] = ['fly0157_rotating_starfield_imaging_T2_trial_1_14602000.abf']
    fly_record['experiments']['img_starfield_t2_rep1']['tiff_file_names'] = ['/T2_trial1/T2_trial1_MMStack.ome.tif']
    fly_record['experiments']['img_starfield_t2_rep1']['sequence_pattern_names'] = starfield_pattern_names_6_0_2014
    fly_record['experiments']['img_starfield_t2_rep1'].create_group('sequences')

    flynum = str(158)
    fly_db.create_group(flynum)
    fly_record =fly_db[flynum]
    fly_record['flynum'] = flynum
    fly_record.create_group('experiments')
    fly_record['experiments'].create_group('img_starfield_t2_rep1')
    #fly_record['experiments'].create_group('b1_azm_expansion_tuning')
    fly_record['experiments']['img_starfield_t2_rep1']['axon_file_names'] = ['fly0158_rotating_starfield_imaging_T2_trial_1_14602002.abf']
    fly_record['experiments']['img_starfield_t2_rep1']['tiff_file_names'] = ['/T2_trial1/T2_trial1_MMStack.ome.tif']
    fly_record['experiments']['img_starfield_t2_rep1']['sequence_pattern_names'] = starfield_pattern_names_6_0_2014
    fly_record['experiments']['img_starfield_t2_rep1'].create_group('sequences')

    flynum = str(159)
    fly_db.create_group(flynum)
    fly_record =fly_db[flynum]
    fly_record['flynum'] = flynum
    fly_record.create_group('experiments')
    fly_record['experiments'].create_group('img_starfield_t2_rep1')
    #fly_record['experiments'].create_group('b1_azm_expansion_tuning')
    fly_record['experiments']['img_starfield_t2_rep1']['axon_file_names'] = ['fly0159_rotating_starfield_imaging_T2_trial_1_14602004.abf']
    fly_record['experiments']['img_starfield_t2_rep1']['tiff_file_names'] = ['/T2_trial1/T2_trial1_MMStack.ome.tif']
    fly_record['experiments']['img_starfield_t2_rep1']['sequence_pattern_names'] = starfield_pattern_names_6_0_2014
    fly_record['experiments']['img_starfield_t2_rep1'].create_group('sequences')

    flynum = str(160)
    fly_db.create_group(flynum)
    fly_record =fly_db[flynum]
    fly_record['flynum'] = flynum
    fly_record.create_group('experiments')
    fly_record['experiments'].create_group('img_starfield_t2_rep1')
    #fly_record['experiments'].create_group('b1_azm_expansion_tuning')
    fly_record['experiments']['img_starfield_t2_rep1']['axon_file_names'] = ['fly0160_rotating_starfield_imaging_T2_trial_1_14602007.abf']
    fly_record['experiments']['img_starfield_t2_rep1']['tiff_file_names'] = ['/T2_trial1/T2_trial1_MMStack.ome.tif']
    fly_record['experiments']['img_starfield_t2_rep1']['sequence_pattern_names'] = starfield_pattern_names_6_0_2014
    fly_record['experiments']['img_starfield_t2_rep1'].create_group('sequences')

    flynum = str(161)
    fly_db.create_group(flynum)
    fly_record =fly_db[flynum]
    fly_record['flynum'] = flynum
    fly_record.create_group('experiments')
    fly_record['experiments'].create_group('img_starfield_t2_rep1')
    #fly_record['experiments'].create_group('b1_azm_expansion_tuning')
    fly_record['experiments']['img_starfield_t2_rep1']['axon_file_names'] = ['fly0161_rotating_starfield_imaging_T2_trial_1_14603000.abf']
    fly_record['experiments']['img_starfield_t2_rep1']['tiff_file_names'] = ['/T2_trial1/T2_trial1_MMStack.ome.tif']
    fly_record['experiments']['img_starfield_t2_rep1']['sequence_pattern_names'] = starfield_pattern_names_6_0_2014
    fly_record['experiments']['img_starfield_t2_rep1'].create_group('sequences')

    flynum = str(162)
    fly_db.create_group(flynum)
    fly_record =fly_db[flynum]
    fly_record['flynum'] = flynum
    fly_record.create_group('experiments')
    fly_record['experiments'].create_group('img_starfield_t2_rep1')
    #fly_record['experiments'].create_group('b1_azm_expansion_tuning')
    fly_record['experiments']['img_starfield_t2_rep1']['axon_file_names'] = ['fly0162_rotating_starfield_imaging_T2_trial_1_14603003.abf']
    fly_record['experiments']['img_starfield_t2_rep1']['tiff_file_names'] = ['/T2_trial1/T2_trial1_MMStack.ome.tif']
    fly_record['experiments']['img_starfield_t2_rep1']['sequence_pattern_names'] = starfield_pattern_names_6_0_2014
    fly_record['experiments']['img_starfield_t2_rep1'].create_group('sequences')

    flynum = str(163)
    fly_db.create_group(flynum)
    fly_record =fly_db[flynum]
    fly_record['flynum'] = flynum
    fly_record.create_group('experiments')
    fly_record['experiments'].create_group('img_starfield_t2_rep1')
    #fly_record['experiments'].create_group('b1_azm_expansion_tuning')
    fly_record['experiments']['img_starfield_t2_rep1']['axon_file_names'] = ['fly0163_rotating_starfield_imaging_T2_trial_1_14603009.abf']
    fly_record['experiments']['img_starfield_t2_rep1']['tiff_file_names'] = ['/T2_trial1/T2_trial1_MMStack.ome.tif']
    fly_record['experiments']['img_starfield_t2_rep1']['sequence_pattern_names'] = starfield_pattern_names_6_0_2014
    fly_record['experiments']['img_starfield_t2_rep1'].create_group('sequences')

    flynum = str(164)
    fly_db.create_group(flynum)
    fly_record =fly_db[flynum]
    fly_record['flynum'] = flynum
    fly_record.create_group('experiments')
    fly_record['experiments'].create_group('img_starfield_t2_rep1')
    #fly_record['experiments'].create_group('b1_azm_expansion_tuning')
    fly_record['experiments']['img_starfield_t2_rep1']['axon_file_names'] = ['fly0164_rotating_starfield_imaging_T2_trial_1_14603012.abf']
    fly_record['experiments']['img_starfield_t2_rep1']['tiff_file_names'] = ['/T2_trial1/T2_trial1_MMStack.ome.tif']
    fly_record['experiments']['img_starfield_t2_rep1']['sequence_pattern_names'] = starfield_pattern_names_6_0_2014
    fly_record['experiments']['img_starfield_t2_rep1'].create_group('sequences')

    flynum = str(165)
    fly_db.create_group(flynum)
    fly_record =fly_db[flynum]
    fly_record['flynum'] = flynum
    fly_record.create_group('experiments')
    fly_record['experiments'].create_group('img_starfield_t2_rep1')
    #fly_record['experiments'].create_group('b1_azm_expansion_tuning')
    fly_record['experiments']['img_starfield_t2_rep1']['axon_file_names'] = ['fly0165_rotating_starfield_imaging_T2_trial_1_14603017.abf']
    fly_record['experiments']['img_starfield_t2_rep1']['tiff_file_names'] = ['/T2_trial1/T2_trial1_MMStack.ome.tif']
    fly_record['experiments']['img_starfield_t2_rep1']['sequence_pattern_names'] = starfield_pattern_names_6_0_2014
    fly_record['experiments']['img_starfield_t2_rep1'].create_group('sequences')

    flynum = str(166)
    fly_db.create_group(flynum)
    fly_record =fly_db[flynum]
    fly_record['flynum'] = flynum
    fly_record.create_group('experiments')
    fly_record['experiments'].create_group('img_starfield_t2_rep1')
    #fly_record['experiments'].create_group('b1_azm_expansion_tuning')
    fly_record['experiments']['img_starfield_t2_rep1']['axon_file_names'] = ['fly0166_rotating_starfield_imaging_T2_trial_1_14605000.abf']
    fly_record['experiments']['img_starfield_t2_rep1']['tiff_file_names'] = ['/T2_trial1/T2_trial1_MMStack.ome.tif']
    fly_record['experiments']['img_starfield_t2_rep1']['sequence_pattern_names'] = starfield_pattern_names_6_0_2014
    fly_record['experiments']['img_starfield_t2_rep1'].create_group('sequences')

    flynum = str(167)
    fly_db.create_group(flynum)
    fly_record =fly_db[flynum]
    fly_record['flynum'] = flynum
    fly_record.create_group('experiments')
    fly_record['experiments'].create_group('img_starfield_t2_rep1')
    #fly_record['experiments'].create_group('b1_azm_expansion_tuning')
    fly_record['experiments']['img_starfield_t2_rep1']['axon_file_names'] = ['fly0167_rotating_starfield_imaging_T2_trial_1_14605001.abf']
    fly_record['experiments']['img_starfield_t2_rep1']['tiff_file_names'] = ['/T2_trial1/T2_trial1_MMStack.ome.tif']
    fly_record['experiments']['img_starfield_t2_rep1']['sequence_pattern_names'] = starfield_pattern_names_6_0_2014
    fly_record['experiments']['img_starfield_t2_rep1'].create_group('sequences')

    flynum = str(168)
    fly_db.create_group(flynum)
    fly_record =fly_db[flynum]
    fly_record['flynum'] = flynum
    fly_record.create_group('experiments')
    fly_record['experiments'].create_group('img_starfield_t2_rep1')
    #fly_record['experiments'].create_group('b1_azm_expansion_tuning')
    fly_record['experiments']['img_starfield_t2_rep1']['axon_file_names'] = ['fly0168_rotating_starfield_imaging_T2_trial_1_14605002.abf']
    fly_record['experiments']['img_starfield_t2_rep1']['tiff_file_names'] = ['/T2_trial1/T2_trial1_MMStack.ome.tif']
    fly_record['experiments']['img_starfield_t2_rep1']['sequence_pattern_names'] = starfield_pattern_names_6_0_2014
    fly_record['experiments']['img_starfield_t2_rep1'].create_group('sequences')

    flynum = str(169)
    fly_db.create_group(flynum)
    fly_record =fly_db[flynum]
    fly_record['flynum'] = flynum
    fly_record.create_group('experiments')
    fly_record['experiments'].create_group('img_starfield_t2_rep1')
    #fly_record['experiments'].create_group('b1_azm_expansion_tuning')
    fly_record['experiments']['img_starfield_t2_rep1']['axon_file_names'] = ['fly0169_rotating_starfield_imaging_T2_trial_1_14605004.abf']
    fly_record['experiments']['img_starfield_t2_rep1']['tiff_file_names'] = ['/T2_trial1/T2_trial1_MMStack.ome.tif']
    fly_record['experiments']['img_starfield_t2_rep1']['sequence_pattern_names'] = starfield_pattern_names_6_0_2014
    fly_record['experiments']['img_starfield_t2_rep1'].create_group('sequences')

    #######
    #######
    #######
    flynum = str(170)
    fly_db.create_group(flynum)
    fly_record =fly_db[flynum]
    fly_record['flynum'] = flynum
    fly_record.create_group('experiments')
    fly_record['experiments'].create_group('img_starfield_t2_rep1')
    #fly_record['experiments'].create_group('b1_azm_expansion_tuning')
    fly_record['experiments']['img_starfield_t2_rep1']['axon_file_names'] = ['fly0170_rotating_starfield_imaging_T2_trial_1_14605004.abf']
    fly_record['experiments']['img_starfield_t2_rep1']['tiff_file_names'] = ['/T2_trial1/T2_trial1_MMStack.ome.tif']
    fly_record['experiments']['img_starfield_t2_rep1']['sequence_pattern_names'] = starfield_pattern_names_6_0_2014
    fly_record['experiments']['img_starfield_t2_rep1'].create_group('sequences')

    flynum = str(171)
    fly_db.create_group(flynum)
    fly_record =fly_db[flynum]
    fly_record['flynum'] = flynum
    fly_record.create_group('experiments')
    fly_record['experiments'].create_group('img_starfield_t2_rep1')
    #fly_record['experiments'].create_group('b1_azm_expansion_tuning')
    fly_record['experiments']['img_starfield_t2_rep1']['axon_file_names'] = ['fly0171_rotating_starfield_imaging_T2_trial_1_14605004.abf']
    fly_record['experiments']['img_starfield_t2_rep1']['tiff_file_names'] = ['/T2_trial1/T2_trial1_MMStack.ome.tif']
    fly_record['experiments']['img_starfield_t2_rep1']['sequence_pattern_names'] = starfield_pattern_names_6_0_2014
    fly_record['experiments']['img_starfield_t2_rep1'].create_group('sequences')

    flynum = str(172)
    fly_db.create_group(flynum)
    fly_record =fly_db[flynum]
    fly_record['flynum'] = flynum
    fly_record.create_group('experiments')
    fly_record['experiments'].create_group('img_starfield_t2_rep1')
    #fly_record['experiments'].create_group('b1_azm_expansion_tuning')
    fly_record['experiments']['img_starfield_t2_rep1']['axon_file_names'] = ['fly0172_rotating_starfield_imaging_T2_trial_1_14605004.abf']
    fly_record['experiments']['img_starfield_t2_rep1']['tiff_file_names'] = ['/T2_trial1/T2_trial1_MMStack.ome.tif']
    fly_record['experiments']['img_starfield_t2_rep1']['sequence_pattern_names'] = starfield_pattern_names_6_0_2014
    fly_record['experiments']['img_starfield_t2_rep1'].create_group('sequences')

    flynum = str(173)
    fly_db.create_group(flynum)
    fly_record =fly_db[flynum]
    fly_record['flynum'] = flynum
    fly_record.create_group('experiments')
    fly_record['experiments'].create_group('img_starfield_t2_rep1')
    #fly_record['experiments'].create_group('b1_azm_expansion_tuning')
    fly_record['experiments']['img_starfield_t2_rep1']['axon_file_names'] = ['fly0173_rotating_starfield_imaging_T2_trial_1_14605004.abf']
    fly_record['experiments']['img_starfield_t2_rep1']['tiff_file_names'] = ['/T2_trial1/T2_trial1_MMStack.ome.tif']
    fly_record['experiments']['img_starfield_t2_rep1']['sequence_pattern_names'] = starfield_pattern_names_6_0_2014
    fly_record['experiments']['img_starfield_t2_rep1'].create_group('sequences')

    flynum = str(174)
    fly_db.create_group(flynum)
    fly_record =fly_db[flynum]
    fly_record['flynum'] = flynum
    fly_record.create_group('experiments')
    fly_record['experiments'].create_group('img_starfield_t2_rep1')
    #fly_record['experiments'].create_group('b1_azm_expansion_tuning')
    fly_record['experiments']['img_starfield_t2_rep1']['axon_file_names'] = ['fly0174_rotating_starfield_imaging_T2_trial_1_14605004.abf']
    fly_record['experiments']['img_starfield_t2_rep1']['tiff_file_names'] = ['/T2_trial1/T2_trial1_MMStack.ome.tif']
    fly_record['experiments']['img_starfield_t2_rep1']['sequence_pattern_names'] = starfield_pattern_names_6_0_2014
    fly_record['experiments']['img_starfield_t2_rep1'].create_group('sequences')

    flynum = str(175)
    fly_db.create_group(flynum)
    fly_record =fly_db[flynum]
    fly_record['flynum'] = flynum
    fly_record.create_group('experiments')
    fly_record['experiments'].create_group('img_starfield_t2_rep1')
    #fly_record['experiments'].create_group('b1_azm_expansion_tuning')
    fly_record['experiments']['img_starfield_t2_rep1']['axon_file_names'] = ['fly0175_rotating_starfield_imaging_T2_trial_1_14605004.abf']
    fly_record['experiments']['img_starfield_t2_rep1']['tiff_file_names'] = ['/T2_trial1/T2_trial1_MMStack.ome.tif']
    fly_record['experiments']['img_starfield_t2_rep1']['sequence_pattern_names'] = starfield_pattern_names_6_0_2014
    fly_record['experiments']['img_starfield_t2_rep1'].create_group('sequences')


    return fly_db

