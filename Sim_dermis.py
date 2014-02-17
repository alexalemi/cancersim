from code import cancer
import pylab as py

from config import config

XSIZE = config['XSIZE']
YSIZE = config['YSIZE']
config['first_cancer_cell_yoffset']=-25
config['pressure_filename']= 'pressure/indaderma_pressure.dat'
config['cancer_evolution_filename']= 'evolution/indaderma_cancer_evolution.dat'

py.figure(1)
py.show()

Q = cancer.CancerSim(config)
Q._setup()

py.ioff()

Q.plot_sized_cells()
py.draw()

imagepath = "2D/indaderma_{:06d}.png"

Q.save('states/indaderma_state_' + str(0) + '.dat')



for i in range(10):
    Q.time_step()
    Q.plot_links()
    py.savefig(imagepath.format(i))
    Q.save('below/states/indaderma_state_' + str(i+1) + '.dat')
    Q.plot_links()
    #py.xlim((5,16))
    #py.ylim((9,17))
    py.draw()
      
    

