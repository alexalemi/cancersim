from code import cancer
from code import cells
import pylab as py

from config import config
import scipy as sp


config['XSIZE'] = 100
config['YSIZE'] = 100
config['basal_height'] = 5
config['first_cancer_cell_yoffset'] = 1.5



def startoff():

    #py.figure(1)
    #py.show()

    Q = cancer.CancerSim(config)

    firstcell = cells.Cell([0,0],Q.cancer,1)

    Q.cancer_cells.append(firstcell)
    Q.add_cell(firstcell)

    #Q.plot_sized_cells()
    #py.draw()

    return Q

#Q._setup()

#py.ioff()

clusterdir = "clusters/"

sizes = [60,80,100]

from uuid import uuid4
import os

def run(steps=40):
    pk = uuid4()
    for i in range(steps):
        Q.time_step()
        if i%10 == 0:
            print i*1./steps
        Q.plot_links()
        py.xlim((-25,25))
        py.ylim((-25,config['YSIZE']))
        py.draw()

    xyfile = os.path.join(clusterdir,str(steps),str(pk) + ".gz")
    #Q.save('clusters/40/{:08d}.dat'.format(pk))
    sp.savetxt(xyfile,Q.get_pos_arr(force=True))
    #py.draw()  


import time
from emailme import send_email

if __name__ == "__main__":

    for size in sizes:
        for r in xrange(1000):
            starttime = time.time()
            Q = startoff()
            run(size)
            elapsed_time = time.time() - starttime
            print "This time took {} seconds".format(elapsed_time)
            send_email(message="The run of size={} finished, it took {} seconds".format(size, elapsed_time))

