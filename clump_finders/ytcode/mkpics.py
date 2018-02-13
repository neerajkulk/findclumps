import glob
import yt
yt.enable_parallelism()
import numpy as np

def mkpic(fname):
    # make a 2048^2 image of simulation domain for making movies (no colorbar or axes)                                                                                            
    ds = yt.load(fname)
    dim = ds.domain_dimensions[0]*(2**ds.max_level)
    slc = yt.SlicePlot(ds, 'z', 'density')
    slc.set_cmap('density', 'binary')
    slc.set_buff_size(2048)
    slc.hide_colorbar()
    slc.hide_axes()
    slc.save(fname + '.png')
    # show a 2048 image of the center of simulation domain                                                                                return

#mkpic("/Volumes/LaCie/simdata/production_runs/thin_8192/shatter.out3.00080.athdf")                                                        


num_procs = 8 # use 64 for stampede                                                                                                        
my_storage = {}

files = glob.glob('/scratch/04325/neerajk/production_runs/32768thick/*athdf')


for sto, file in yt.parallel_objects(files, num_procs, storage = my_storage):
    mkpic(file)




    



