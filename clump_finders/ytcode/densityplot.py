import glob
import yt
yt.enable_parallelism()
import numpy as np

def densplot(fname):
    # make a 1024^2 image of simulation domain                                                                                             
    ds = yt.load(fname)
    dim = ds.domain_dimensions[0]*(2**ds.max_level)
    slc = yt.SlicePlot(ds, 'z', 'density')
    slc.set_cmap('density', 'binary')
    slc.set_buff_size(1024)
    slc.save(fname + '.png')
    # show a 2048 image of the center of simulation domain                                                                                 
    slc.set_buff_size(2048)
    slc.zoom(dim/2048)
    slc.save(fname + '.zoom.png')
    return

#mkpic("/Volumes/LaCie/simdata/production_runs/thin_8192/shatter.out3.00080.athdf")                                                        


num_procs = 8 # use 64 for stampede                                                                                                        
my_storage = {}

files = glob.glob('/scratch/04325/neerajk/production_runs/32768thick/*athdf')


for sto, file in yt.parallel_objects(files, num_procs, storage = my_storage):
    densplot(file)




    



