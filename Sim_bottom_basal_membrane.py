from code import cancer
import pylab as py

from config import config

XSIZE = config['XSIZE']
YSIZE = config['YSIZE']
config['first_cancer_cell_xoffset']=-1.0*config['XSIZE']/(3.0*config['basal_wavenumber'])
config['pressure_filename']= 'pressure/bottom_basal_membrane_pressure.dat'
config['cancer_evolution_filename']= 'evolution/bottom_basal_membrane_cancer_evolution.dat'

py.figure(1)
py.show()

Q = cancer.CancerSim(config)
Q._setup()

py.ioff()

Q.plot_sized_cells()
py.draw()

imagepath = "2D/bottom_basal_membrane_{:06d}.png"

Q.save('states/bottom_basal_membrane_state_' + str(0) + '.dat')

for i in range(10):
    Q.time_step()
    Q.plot_links()
    py.savefig(imagepath.format(i))
    Q.save('states/bottom_basal_membrane_state_' + str(i+1) + '.dat')
    Q.plot_links()
    #py.xlim((5,17))
    #py.ylim((9,17))
    py.draw()
      
    

