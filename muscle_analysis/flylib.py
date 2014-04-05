import numpy as np

def get_signal_dict(filename):
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
