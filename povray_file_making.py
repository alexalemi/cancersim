import numpy as np
#import numarray as numar
from numpy.random import *
from scipy import *
from scipy.io import *
from scipy.spatial.distance import euclidean
from numpy import *
import matplotlib  
import pylab  
import os                         # For issuing commands to the OS.
import sys                        # For determining the Python version.
import cPickle as pickle 
import pylab as py
import matplotlib.cm as cm
import matplotlib.colors as colors
import numexpr as ne
from code import cancer
from helper import norm
from config import config




def povray_making_movie(position) :

    X_size = 60

    # position == 0 ----> dermis
    # position == 1 ----> epidermis
    # position == 2 ----> top of basal membrane
    # position == 3 ----> bottom of basal membrane
    #position = 1

   
    if position == 0:
        open_file = 'states/indaderma_state_'
        open_file_pressure = 'pressure/indaderma_pressure.dat'
        write_file = '3D/indaderma_3D_'
        
    if position == 1:
        open_file = 'states/indaepidermis_state_'
        open_file_pressure = 'pressure/indaepidermis_pressure.dat'
        write_file = '3D/indaepidermis_3D_'
       
    if position == 2:
        open_file = 'states/top_basal_membrane_state_'
        open_file_pressure = 'pressure/top_basal_membrane_pressure.dat'
        write_file = '3D/top_basal_membrane_3D_'
        
    if position == 3:
        open_file = 'states/bottom_basal_membrane_state_'
        open_file_pressure = 'pressure/bottom_basal_membrane_pressure.dat'
        write_file = '3D/bottom_basal_membrane_3D_'
    


    filename_upload = open_file_pressure
    upload  = loadtxt(filename_upload)
    stress = empty((len(upload), 1))
    for i in range(0,len(upload)) :
        stress[i] = upload[i][1]
    
    max_number_of_cells = int(1 + upload[len(upload)-1][0])
    print 'max number of cells = ', max_number_of_cells

    for num_of_cell in range (0, max_number_of_cells):
        print 'num_of_cell = ', num_of_cell
        open_file_state = open_file + str(num_of_cell) +'.dat'
        with open(open_file_state,'r') as f:
            config, cells, links, ghosts, T = pickle.load(f)
        stress_partial = empty((num_of_cell + 1, 1))
        for i in range (0, num_of_cell+1) :
            stress_partial[i] = stress[(num_of_cell+1)*num_of_cell/2+i]
        print stress_partial
        if len(stress_partial)>1 :
            stress_partial = (stress_partial-stress_partial.min())/(stress_partial.max()-stress_partial.min())*0.9+0.1
        else :
            stress_partial[0] = 0.5
        col_array = []
        for i in range(0,len(stress_partial)) :
            rgb_color =  cm.hot(1-stress_partial[i],1.0)
            col_array.append(rgb_color)

        write_file_3D =  write_file +  str(num_of_cell) + '.txt'
        file_povray=open(write_file_3D,'w')
 
        numero_cancer = 0
        for cell in cells:
            if cell.type.name == 'tDermal':
                color = 'color LimeGreen'
            if cell.type.name == 'Epidermal':
                color = 'color MediumBlue'
            if cell.type.name == 'Basal':
                color = 'color Gray20'
            if cell.type.name == 'Corneum':
                color = 'color MediumVioletRed'
            if cell.type.name == 'Cancer':
                color = 'color <' + str(col_array[numero_cancer][0][0]) + ',' + str(col_array[numero_cancer][0][1]) + ',' +  str(col_array[numero_cancer][0][2]) + '>'
                numero_cancer = numero_cancer + 1
            s = 'sphere { <0, 0, 0>,' + str(cell.radius) + ' material {texture {pigment {' + color + '} finish { specular specularvalue roughness roughnessvalue reflection phongreflection }}}translate<' + str(cell.pos[0]) + ',' + str(cell.pos[1]) + ', 0.0 >}\n'
            file_povray.write(s)

        cutoff = X_size/2
        for link in links:
            if (link.two.pos[0]!=link.one.pos[0]) and (link.two.pos[1]!=link.one.pos[1]):
                d12 = link.one.pos-link.two.pos
                abs_d12 = norm(d12)
                if abs_d12 < cutoff :
                    color = 'color White'
                    s = 'cylinder { <' + str(link.one.pos[0]) + ',' + str(link.one.pos[1]) + ', 0.0 >,<' + str(link.two.pos[0]) + ',' + str(link.two.pos[1]) + ', 0.0 >,'  + str(0.1) + ' material {texture {pigment {' + color + '} finish { phong phongvalue_cyl phong_size phongsize_cyl reflection phongreflection_cyl}}}}\n'
                    file_povray.write(s)
        file_povray.close

    
    
    
    
      
    
