import yt
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
    my_grid = ds.covering_grid(level=ds.max_level,
                               left_edge=ds.domain_left_edge,
                               dims=ds.domain_dimensions,
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
    clumpdatat = [ sliceclumps(slice) for slice in np.transpose(tdata) ]  #map slice clumps over thresholded data
    
#    allclumps = np.concatenate((clumpdatat,clumpdata))


    flat_listx = [item for sublist in clumpdata for item in sublist] #flatten the list
    flat_listy = [item for sublist in clumpdatat for item in sublist] #flatten the list

    cx = Counter(flat_listx)
    cy = Counter(flat_listy)

    # write clumps to file
    filenamex = 'outputx.clumps'
    filenamey = 'outputy.clumps'
    
    with open (filenamex, 'w') as f:
        f.write("# time = {0}\n".format(c_time))
        f.write("# res = {0}\n".format(resolution))
        f.write("# mintemp = {0}\n".format(tempdata.min()))
        f.write("# [1] = size, [2] = count\n")
        for (size, count) in cx.items():
            f.write("{0}\t{1}\n".format(size, count))

    with open (filenamey, 'w') as f:
        f.write("# time = {0}\n".format(c_time))
        f.write("# res = {0}\n".format(resolution))
        f.write("# mintemp = {0}\n".format(tempdata.min()))
        f.write("# [1] = size, [2] = count\n")
        for (size, count) in cy.items():
            f.write("{0}\t{1}\n".format(size, count))


    return




findclumps("/Volumes/LaCie/simdata/test_run/shatter.out3.00028.athdf")
            
            







# BUNCH OF ARCHIVED CODE!!!!!!

   # flat_list = [item for sublist in tdata for item in sublist]
   # print('number of zeros and ones')
   # print(Counter(flat_list))

    # clumpdata = [ sliceclumps(slice) for slice in tdata ]  #map slice clumps over thresholded data
    # print(length(clumpdata))
    # clumpdatat = [ sliceclumps(slice) for slice in np.transpose(tdata) ]  #map slice clumps over thresholded data
    # allclumps = np.concatenate((clumpdatat,clumpdata))
    # flat_list = [item for sublist in allclumps for item in sublist] #flatten the list
    # c = Counter(flat_list)
    # filename = 'output.clumps'
    
    # with open (filename, 'w') as f:
    #     f.write("# time = {0}\n".format(c_time))
    #     f.write("# res = {0}\n".format(resolution))
    #     f.write("# [1] = size, [2] = count\n")
    #     for (size, count) in c.items():
    #         f.write("{0}\t{1}\n".format(size, count))



    
# NOW START SLICING
#data = map(sliceclumps,t)
#data = [ sliceclumps(slice) for slice in t ]  #map slice clumps over thresholded data
#flat_list = [item for sublist in data for item in sublist] #flatten the list

# from collections import Counter
# c = Counter(flat_list)

# # read these from file
# time = 4.0
# resolution = 4096

# filename = 'output.clumps'

# with open (filename, 'w') as f:
#     f.write("# time = {0}\n".format(time))
#     f.write("# res = {0}\n".format(resolution))
#     f.write("# [1] = size, [2] = count\n")
#     for (size, count) in c.items():
#         f.write("{0}\t{1}\n".format(size, count))


    # t = temperature
    # tmin = 0.1
    # t[t < tmin] = 1.0
    # t[t > tmin] = 0.0
    # num_zeros = (t == 0.0).sum()
    # num_ones = (t == 1.0).sum()
    # t = t[0:,0:,0]
    # print(t.shape)
    # np.savetxt("%s_density" % ds.basename, density)
    # t=np.array([0,0,1,1,1,0,0,0,1,1,1,1,1,1,0,0])
    # lw, num = measurements.label(t)
    # print(lw)
    # area = measurements.sum(t, lw, index=range(lw.max() + 1))
    # print(area)
    # print(t)

    # look for scikit-image resize with antialiasing

    
import glob
import os

files = glob.glob("*.athdf")

for file in files:
    pre, ext = os.path.splittext(file)
    if not (file_exists(pre + '.clumps')):
        analyze(file)
