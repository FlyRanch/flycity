import neo
import quantities as pq
import h5py
import flylib
import db_access as dba
import numpy as np
i = 17
imaging_group = [151,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169]
fly_db = dba.get_db()
swarm = flylib.Squadron(fly_db,imaging_group)
expmnt = swarm.flies[i].experiments[u'img_starfield_t2_rep1']
print(imaging_group[i])
import imaging_viewer
viewer = imaging_viewer.ImagingViewer(expmnt)
viewer.run()

expmnt.set_roi_data(viewer.muscle_data)
fly_db.flush()
fly_db.close()