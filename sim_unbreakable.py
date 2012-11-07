from code import cancer
import pylab as py

from config import config

XSIZE = config['XSIZE']
YSIZE = config['YSIZE']
config['first_cancer_cell_yoffset']=1
config['basal_cell_params']['maxstretch'] = 10000


py.figure(1)
#py.show()

Q = cancer.CancerSim(config)
Q._setup()

py.ioff()

Q.plot_sized_cells()
py.draw()


basepath = "pics/unbreakable2/"
imagepath = basepath + "{:06d}.png"


for i in range(1000):
    Q.time_step()
    #Q.plot_links()
    Q.save(basepath + 'state.dat')
    #Q.plot_links()
    py.xlim((-5,XSIZE+5))
    py.ylim((-5,YSIZE+5))
    py.savefig(imagepath.format(i))
    #py.draw()
      
    

