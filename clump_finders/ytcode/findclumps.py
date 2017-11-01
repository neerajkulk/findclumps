import glob
import os.path
import yt
yt.enable_parallelism()
import numpy as np
from scipy.ndimage import measurements

from collections import Counter

#define simulation parameters
tlim = 15.0
alpha = -1.0
gm = 1.66666666666666666666666666667
drat = 10.0

f = gm**1.5/((gm - 1.0)*(2.0-alpha))

lambda0 = 1/f
lambda0 *= 1/(drat**(2.5-alpha))
lambda0 *= 4.0


#define relavant functions

def cstcool(x):
    times = [2**(-i-1.0) for i in range(5)]
    times[:] = [i / sum(times) for i in times]
    times = np.cumsum(times)
    times=np.append(times,x)
    times.sort()
    loc = np.where( times==x )
    loc = loc[0].flat[0]
    return 4.0**(-loc-1.0)

def tfloor(x):
    return (f*lambda0*cstcool(x))**(2.0/(5.0-2.0*alpha))

def get_time(filename):
    ds = yt.load(filename)
    return ds.current_time



def temp_plus(filename):
    ds = yt.load(filename)
    amr_factor = np.array([2**ds.max_level,2**ds.max_level,1])
    my_grid = ds.covering_grid(level=ds.max_level,
                               left_edge=ds.domain_left_edge,
                               dims=ds.domain_dimensions * amr_factor,
                               fields=["density", "pressure"])
    temp=my_grid["pressure"]/my_grid["density"]
    return np.array(temp)[0:,0:,0]



def mybinarize(data,tflr):
    tdata = data < 2.0*tflr
    tdata = tdata.astype("int")
    return tdata


def sliceclumps(slice):
    lw, num = measurements.label(slice)
    area = measurements.sum(slice, lw, index=range(lw.max() + 1))
    area = np.delete(area,0)
    return area.astype(np.int)


# Main function to analyze files
def findclumps(fname):
    tempdata = temp_plus(fname)
    c_time = get_time(fname)
    tdata = mybinarize(tempdata, tfloor(c_time/tlim))
    resolution = tdata.shape[0]

    clumpdata = [ sliceclumps(slice) for slice in tdata ]  #map slice clumps over thresholded data
    clumpdata_t = [ sliceclumps(slice) for slice in np.transpose(tdata) ]  #map slice clumps over thresholded data
    
    flat_list = [item for sublist in clumpdata for item in sublist] #flatten the list
    flat_list_t = [item for sublist in clumpdata_t for item in sublist] #flatten the list

    c = Counter(flat_list)
    c_t = Counter(flat_list_t)


    # write clumps to file                                                                                           
    filename = fname + '.clumps'
    filename_t = fname + '.t.clumps'

    with open (filename, 'w') as f:
        f.write("# time = {0}\n".format(c_time))
        f.write("# res = {0}\n".format(resolution))
        f.write("# mintemp = {0}\n".format(tempdata.min()))
        f.write("# [1] = size, [2] = count\n")
        for (size, count) in c.items():
            f.write("{0}\t{1}\n".format(size, count))
        f.write("all_clumps_found")
        f.close()

    with open (filename_t, 'w') as f:
        f.write("# time = {0}\n".format(c_time))
        f.write("# res = {0}\n".format(resolution))
        f.write("# mintemp = {0}\n".format(tempdata.min()))
        f.write("# [1] = size, [2] = count\n")
        for (size, count) in c_t.items():
            f.write("{0}\t{1}\n".format(size, count))
        f.write("all_clumps_found")
        f.close()
    return




#findclumps('/Volumes/LaCie/simdata/test_run/shatter.out3.00028.athdf')                                              

#parallelize and map over all clumps                                                                                 

num_procs = 64 # use 64 for stampede                                                                                 
my_storage = {}

files = glob.glob('/scratch/04325/neerajk/production_runs/thin_8192/*.athdf')


for sto, file in yt.parallel_objects(files, num_procs, storage = my_storage):
    if not (os.path.isfile(file + '.clumps') and os.path.isfile(file + '.t.clumps')):
        findclumps(file)
