import glob
import os.path
import yt
yt.enable_parallelism()
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
import yt


def temp_plus(filename):
    ds = yt.load(filename)
    amr_factor = np.array([2**ds.max_level,2**ds.max_level,1])
    my_grid = ds.covering_grid(level=ds.max_level,
                               left_edge=ds.domain_left_edge,
                               dims=ds.domain_dimensions * amr_factor,
                               fields=["density", "pressure"])
    temp=my_grid["pressure"]/my_grid["density"]
    return np.array(temp)[0:,0:,0]



def write_hist(fname):
    data = temp_plus(fname)
    flatdata = np.ndarray.flatten(data)
    c = Counter(np.round(flatdata, 2))

    with open(fname + '.tmp', 'w') as f:
        f.write('# [1] = temp, [2] = PDF (unnormalized)\n')
        for x, y in sorted(c.items()):
            f.write("{0:6.2f}\t{1:6d}\n".format(x, y))
        f.write("all_done")
        f.close()
    return


#write_hist("/Volumes/LaCie/simdata/production_runs/8192thick/shatter.out3.00030.athdf")



#parallelize and map over all files                                                                                 

num_procs = 16
my_storage = {}

files = glob.glob('/scratch/04325/neerajk/production_runs/thin_8192/*.athdf')


for sto, file in yt.parallel_objects(files, num_procs, storage = my_storage):
    if not (os.path.isfile(file + '.tmp')):
        write_hist(file)

