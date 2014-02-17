#Cancer Sim

from numpy import *

import scipy as sp
import pylab as py
import math

import matplotlib.cm as cm
import matplotlib.colors as colors

import cPickle as pickle

from scipy.spatial.distance import euclidean
from math import pow
from scipy.spatial import Delaunay
#from scipy.spatial import KDTree 
from scipy.spatial import cKDTree
from hybridKDTree import KDTree
import random
import time
import pprint

#XSIZE = 20
#YSIZE = 20

from neighborlist import NeighborList

from helper import norm, unitize, disp_func, unitize_arr
import links, cells

from logger import logger
base_logger = logger.getChild('cancer')
base_logger.info('Inside the cancer.py module')

########################################################
### Simulation Class ###################################
########################################################

#try to speed things up a little bit
from scipy import zeros_like, nan_to_num, allclose

import numexpr as ne


import os



if 'CANCERC' in os.environ:
    CANCERC = True
    #import pyximport
    #pyximport.install()
    from forcefunccelltypes import force_func_hertz, force_func_basal, norm, disp_func
    base_logger.info('CYTHON SUPPORT')
else:
    CANCERC = False

    force_func_basal = None
    force_func_hertz = None

    base_logger.info('NO CYTHON SUPPORT')




class CancerSim:
    """ 
        The main Cancer Simulation Class.
        
        Creates an array of Cells, allows for the designation of cancer cells
        And the evolution of the cells thereafter.
    """
    def __init__(self,config):
        """ Initialize the simulation """


        #load the configs
        self.config = config
        self.XSIZE = config['XSIZE']
        self.YSIZE = config['YSIZE']
        self.boxsize = (self.XSIZE,self.YSIZE)

        if config['seed'] is None:
            self.seed = int(time.time())
        else:
            self.seed = config['seed']
            

        self.xi = config['force_cutoff']
        self.a = config['force_magnitude']
        self.basalstrength = config['force_magnitude_basal']
        self.basalcutoff = config['force_cutoff_basal']
        self.basal_height = config['basal_height']
        self.basal_wavenumber = config['basal_wavenumber']
        self.basal_amplitude = config['basal_amplitude']
        self.pressure_filename = config['pressure_filename']
        self.cancer_evolution_filename = config['cancer_evolution_filename']

        sp.random.seed(self.seed)
        random.seed(self.seed)


        #KDTree
        #self._kdtree = None
        #self._kdtree_cache_T = -1


        self._updated = True
        self.T = 0

        # cell types (should be arguments)
        self.cancer = cells.CellType(**config['cancer_cell_params'])
        self.epidermal = cells.CellType(**config['epidermal_cell_params'])
        self.basal = cells.CellType(**config['basal_cell_params'])
        self.dermal = cells.CellType(**config['dermal_cell_params'])
        self.corneum = cells.CellType(**config['stratum_corneum_cell_params'])

        self.num_cells = 0

        # containers
        self.links = links.Links()
        self._cell_arr = sp.array([])

        self.cells = []
        self._ghosts = []
        self._ghost_cutoff = 4
        self._ghost_offset = sp.array([self.boxsize[0],0.])
        self.cancer_cells = []

        self.logger = base_logger.getChild('CancerSim')
        self.logger.info('Initializing CancerSim')

        self.neighs = None

    def _setup(self):

        self._triang_lattice()
        self.jiggle(sigma=self.config['jiggle_sigma'])
        self.delaunay()
        self._freeze_links()
        XSIZE, YSIZE = self.boxsize

        period = 2*3.141592*self.basal_wavenumber/XSIZE
        

        self.add_cancer_cell([self.XSIZE/2.+self.config['first_cancer_cell_xoffset'], self.basal_height + self.basal_amplitude*sin((self.XSIZE/2+self.config['first_cancer_cell_xoffset'])*period) + self.config['first_cancer_cell_yoffset']], self.config['first_cancer_cell_radius'])
        

    def _triang_lattice(self):
        """ Create a triangular grid of points """
        XSIZE, YSIZE = self.boxsize
        period = 2*3.141592*self.basal_wavenumber/XSIZE

        self.logger.info('Setting up the Triangular Lattice...')
        #setup the epicells
        epispacing = self.epidermal.L
        xspace,yspace = epispacing , epispacing * sp.sqrt(3)
        for i in sp.arange(0,XSIZE,xspace):
            for ind,j in enumerate(sp.arange(self.basal_height-self.basal_amplitude+5.0*self.basalcutoff,YSIZE,yspace)):
                if ind:
                    pass
                if j >= self.basal_height+self.basal_amplitude*sin(i*period)+5.0*self.basalcutoff :
                    cell1 = cells.Cell([i,j],self.epidermal,self.num_cells)
                #print 'added epicell at', i, j
                    self.add_cell(cell1)
                if (j+0.5*yspace) > self.basal_height+self.basal_amplitude*sin((i+0.5*xspace)*period) :
                    cell2 = cells.Cell([i+0.5*xspace,j+0.5*yspace],self.epidermal,self.num_cells)
                #print 'added epicell at', i+0.5*xspace, j+0.5*yspace
                    self.add_cell(cell2)

                #add ghosts for first few layers
                if i<self._ghost_cutoff:
                    if ind:
                        if j >= self.basal_height+self.basal_amplitude*sin(i*period)+5.0*self.basalcutoff :
                            ghost1 = cells.GhostCell(cell1,XSIZE,1)
                            self._ghosts.append(ghost1)
                    if (j+0.5*yspace) > self.basal_height+self.basal_amplitude*sin((i+0.5*xspace)*period) :
                        ghost2 = cells.GhostCell(cell2,XSIZE,1)
                        self._ghosts.append(ghost2)

                #add ghosts for last few layers
                if i>(XSIZE-self._ghost_cutoff):
                    if ind:
                        if j >= self.basal_height+self.basal_amplitude*sin(i*period)+5.0*self.basalcutoff :
                            ghost1 = cells.GhostCell(cell1,XSIZE,-1)
                            self._ghosts.append(ghost1)
                    if (j+0.5*yspace) > self.basal_height+self.basal_amplitude*sin((i+0.5*xspace)*period) :
                        ghost2 = cells.GhostCell(cell2,XSIZE,-1)
                        self._ghosts.append(ghost2)

        #setup the bottom cells
        dermalspacing = self.dermal.L
        xspace,yspace = dermalspacing , dermalspacing*sp.sqrt(3)
        for i in sp.arange(0,XSIZE,xspace):
            for ind,j in enumerate(sp.arange(self.basal_height+self.basal_amplitude-5.0*self.basalcutoff,0,-yspace)):
                if j<= self.basal_height+self.basal_amplitude*sin(i*period)-5.0*self.basalcutoff :
                    cell1 = cells.Cell([i,j],self.dermal,self.num_cells)
                #print 'added dermacell at', i, j
                    self.add_cell(cell1)
                if ind and (j+0.5*yspace) <= self.basal_height+self.basal_amplitude*sin((i+0.5*xspace)*period)-5.0*self.basalcutoff:
                    cell2 = cells.Cell([i+0.5*xspace,j+0.5*yspace],self.dermal,self.num_cells)
                    #print 'added dermacell at', i+0.5*xspace, j+0.5*yspace
                    self.add_cell(cell2)

                #add ghosts for first few layers
                if i<self._ghost_cutoff:
                    if j<= self.basal_height+self.basal_amplitude*sin(i*period)-5*self.basalcutoff :
                        ghost1 = cells.GhostCell(cell1,XSIZE,1)
                        ghost2 = cells.GhostCell(cell2,XSIZE,1)
                        self._ghosts.extend([ghost1,ghost2])

                #add ghosts for last few layers
                if i>(XSIZE-self._ghost_cutoff):
                    if j<= self.basal_height+self.basal_amplitude*sin(i*period)-5.0*self.basalcutoff :
                        ghost1 = cells.GhostCell(cell1,XSIZE,-1)
                        ghost2 = cells.GhostCell(cell2,XSIZE,-1)
                        self._ghosts.extend([ghost1,ghost2])

        #setup the middle cells
        basalspacing = self.basal.L
        for i in sp.arange(0,XSIZE,basalspacing/2):
            cell = cells.Cell([i,self.basal_height+self.basal_amplitude*sin(i*period)],self.basal,self.num_cells)
            #print 'added basalcell at', i, self.basal_height+self.basal_amplitude*sin(i*period)
            self.add_cell(cell)
            if i<self._ghost_cutoff:
                ghost = cells.GhostCell(cell,XSIZE,1)
                self._ghosts.append(ghost)
            if i>(XSIZE-self._ghost_cutoff):
                ghost = cells.GhostCell(cell,XSIZE,-1)
                self._ghosts.append(ghost)

        #setup the corneum cells
        corneumspacing = self.corneum.L
        for i in sp.arange(0,XSIZE,corneumspacing):
            cell = cells.Cell([i,YSIZE+2.0*self.basalcutoff],self.corneum,self.num_cells)
            #print 'added corneumcell at', i, YSIZE
            self.add_cell(cell)
            if i<self._ghost_cutoff:
                ghost = cells.GhostCell(cell,XSIZE,1)
                self._ghosts.append(ghost)
            if i>(XSIZE-self._ghost_cutoff):
                ghost = cells.GhostCell(cell,XSIZE,-1)
                self._ghosts.append(ghost)




        self.logger.info('Set up the Triangular Lattice')

    
    def get_pos_arr(self,force=False):
        """ Get an array of all of the cell positions """
        #if self._updated is False or force:
        #    return self._cell_arr
        
        self._cell_arr = sp.zeros((len(self.cells),2))
        for (i,cell) in enumerate(self.cells):
            self._cell_arr[i] = cell.pos
        
        self._updated = False
        return self._cell_arr

    def get_radius_arr(self):
        rad_arr=sp.zeros(len(self.cells))
        for (i,cell) in enumerate(self.cells):
            rad_arr[i] = cell.radius
        return rad_arr

    def _get_kdtree(self,force=False,new=True):
        """ Generate a KDTree for the cells, 
            allows for efficient geometric neighbor computation """
        #if new or self._kdtree_cache_T != self.T or self._updated:
        pos = self.get_pos_arr(force).copy()
        _kdtree = KDTree(pos)

        return _kdtree

    def _get_ckdtree(self,force=False):
        """ Generate a cKDTree """
        pos = self.get_pos_arr(force).copy()
        return cKDTree(pos)

    def _query_point(self,x,r,eps=None):
        """ Get all of the cell inds near point, with radius r """        
        kdtree = self._get_kdtree()

        if eps:
            cell_inds = kdtree.query_ball_point(x,r,eps)
        else:
            cell_inds = kdtree.query_ball_point(x,r)
        cells = [ self.cells[ind] for ind in cell_inds ]
        return cells

    def _get_vel_arr(self):
        """ Get an array of all of the cell velocities """
        vel_arr = sp.zeros((self.num_cells,2))
        for (i,cell) in enumerate(self.cells):
            vel_arr[i] = cell.vel
        return vel_arr

    def _update_pos(self,pos_arr):
        """ Update all of the cell positions with an array """
        for (pos,cell) in zip(pos_arr,self.cells):
            #enact the periodic boundary conditions
            pos[0] = pos[0]%self.XSIZE
            cell.pos = pos
        self._cell_arr = pos_arr
        #self._updated = True

    def _update_vel(self,vel_arr):
        """ Update all of the cell velocities with an array """
        for (vel,cell) in zip(vel_arr,self.cells):
            cell.vel = vel
        
    def _get_ghost_pos_arr(self):
        """ Get all of the ghost positions """
        arr = sp.zeros((len(self._ghosts),2))
        for ind,cell in enumerate(self._ghosts):
            arr[ind] = cell.pos
        return arr

    def _update_ghosts(self):
        """ Update the positions of all of the ghost cells """
        for ghost in self._ghosts:
            ghost.update()

    def jiggle(self,sigma=0.1,ghosts=True):
        """ Jiggle the atom positions """
        pos = self.get_pos_arr()

        sigarr = sp.array([cell.type.L for cell in self.cells])
        randn = sp.randn(self.num_cells,2)

        newpos = pos + sigma*(sigarr*randn.T).T

        self._update_pos(newpos)
        self._updated = True

        if ghosts:
            self._update_ghosts()

        self.logger.info('Jiggled the atoms')

    def _set_radii(self):
        """ set radii as the average of the links starting from each cell """
        for cell in [cell for cell in self.cells if cell.type == self.epidermal]:
            average_length=0.0
            count=0.
            for neigh in self.links.get_neighbors(cell):
                average_length += self.links.get_link(cell,neigh).L/2.0
                count += 1.

            if count:
                cell.radius=average_length/count
            
        for cell in [cell for cell in self.cells if cell.type == self.dermal]:
            cell.radius=self.epidermal.L/2.0

    def _set_radii_min(self):
        """ set radii as the smallest link size """
        for cell in [cell for cell in self.cells if cell.type == self.epidermal]:
            min_length = min([link.L/2. for link in self.links.get_links(cell)])
            #rint min_length

            cell.radius=min_length
            
        for cell in [cell for cell in self.cells if cell.type == self.dermal]:
            cell.radius=self.epidermal.L/2.0

    def _freeze_links(self):
        """ Adjust all of the links to be their current extension """
        for link in self.links:
            link.L = link.extension_without_breaking()
            if (link.one.type.name == 'Dermal'):
                 if (link.two.type.name == 'Basal') :
                     print link.one, link.two, link.L

            if (link.one.type.name == 'Epidermal'):
                 if (link.two.type.name == 'Basal') :
                     print link.one, link.two, link.L

        self._set_radii_min()

        self.logger.info('Froze the links in place')


    def _filter_ghosts(self,one,two):
        if isinstance(one,cells.GhostCell) and isinstance(two,cells.GhostCell):
            raise Exception("DoubleGhost")
        elif isinstance(one,cells.GhostCell):
            return one.original,two
        elif isinstance(two,cells.GhostCell):
            return one,two.original
        else:
            return one,two


    def _clear_links(self):
        """ Clear all Links """
        self.links = links.Links()

    def delaunay(self):
        """ Delaunay routine, sets the initial links """

        self.logger.debug('Running the Delaunay routine')

        #first get the positions of all the cells and the ghosts
        num_cells = len(self.cells)
        num_ghosts = len(self._ghosts)
        fulllist = self.cells + self._ghosts
        num_full = len(fulllist)

        arr = sp.zeros((num_full,2))
        for ind,cell in enumerate(fulllist):
            arr[ind] = cell.pos
        
        #get the Delaunay construction
        tri = Delaunay(arr)

        #add the links
        for i,j,k in tri.vertices:
            cellone = fulllist[i]
            celltwo = fulllist[j]
            cellthree = fulllist[k]
            length_of_bond = norm(cellone.pos - celltwo.pos)
            expected_length = 0.5*(cellone.type.L + celltwo.type.L)
            if length_of_bond < 2*expected_length:
                try:
                    one,two = self._filter_ghosts(cellone,celltwo)
                    self.add_bond(one,two)
                except Exception, e:
                    if e.message=="DoubleGhost":
                        pass
                    else:
                        raise
                try:
                    one,two = self._filter_ghosts(celltwo,cellthree)
                    self.add_bond(one,two)
                except Exception, e:
                    if e.message=="DoubleGhost":
                        pass
                    else:
                        raise
                try:
                    one,two = self._filter_ghosts(cellthree,cellone)
                    self.add_bond(one,two)
                except Exception, e:
                    if e.message=="DoubleGhost":
                        pass
                    else:
                        raise
        
    def add_cell(self,cell):
        """ Add the cell: cell """
        self.cells.append(cell)
        self.num_cells += 1
        self._updated = True
        self.logger.debug('Adding the cell {cell}'.format(cell=cell))

    def add_bond(self,one,two):
        """ Add a bond between cells one and two """
        self.links.add_link(one,two,xsize=self.XSIZE)

        self.logger.debug('Adding a bond between {one} and {two}'.format(one=one,two=two))

    def remove_bond(self,one,two):
        """ Remove a bond between cells one and two """
        self.links.remove_link(one,two)

        self.logger.debug('Removed the link between {one} and {two}'.format(one=one,two=two))

    def remove_cell(self,cell):
        """ Remove the cell: cell, and all bonds for that cell """
        self.cells.remove(cell)
        self.links.remove_cell(cell)

        self.logger.debug('Removed the cell {cell}'.format(cell=cell))

    def get_neighbors(self,cell):
        """ Get the linked neighbor cells of cell """
        return self.links.get_neighbors(cell)

    def add_cancer_cell(self,x,r,eps=None):
        file=open(self.cancer_evolution_filename,'a')
        """ randomly make a cell a cancer cell """
        cells = self._query_point(x,r,eps)
        cells = [cell for cell in cells if cell.type != self.basal]
        if cells:
            cell = random.choice(cells)
            self.cancer_cells.append(cell)
            self.links.remove_cell(cell)

            cell.type = self.cancer

            s =  str(cell.pos[0]) + ' ' + str(cell.pos[1])  + '\n'
            file.write(s)

            self.logger.info('Added a cancer cell: {cell}'.format(cell=cell))
            self._updated = True
        else:
            raise Exception("No targets found at {} within radius {}".format(x,r))

        file.close

    def duplicate_cancer_cell(self,cancer=None,disp_frac = 0.01):
        """ Duplicate the cancer cell: cancer """
        if cancer is None:
            cancer = random.choice(self.cancer_cells)
        
        file=open(self.cancer_evolution_filename,'a')

        self.logger.info('Duplicating a cancer cell...')
        #need to choose a random direction and do the relaxation
        L = disp_frac * cancer.type.L
        theta = sp.rand()*2*sp.pi

        disp = L * sp.array([sp.sin(theta),sp.cos(theta)])

        newcell = cells.Cell(cancer.pos + disp,self.cancer,self.num_cells)
        newcell.radius = cancer.radius

        cancer.pos = cancer.pos - disp

        s =  str(cancer.pos[0]) + ' ' + str(cancer.pos[1])  + '\n'
        file.write(s)
        
        self.cancer_cells.append(newcell)
        self.add_cell(newcell)

        """
        neighs = self.links.get_neighbors(cancer).copy()

        
        for neigh in neighs:
            link_disp =  neigh.pos - cancer.pos
            if sp.vdot(link_disp,disp) >= 0:
                #remove old link, create new one.
                self.links.remove_link(cancer,neigh)
                self.links.add_link(newcell,neigh)
        """

        #self.links.add_link(newcell,cancer)

        self._updated = True
        file.close










    def time_step(self):
        """ Run a time step, duplicate a cancer cell,
             do a FIRE relaxation, and plot """
        
        self.logger.info('Running a time step')
        self.duplicate_cancer_cell()
        self.fire()
        self.plot_sized_cells()
        self.T += 1



    def plot_cells(self,clf=True,fignum=1,ghosts=False,*args,**kwargs):
        """ Plot the current configuration """

        self.logger.info('Plotting the cells')
        pos_arr = self.get_pos_arr()

        py.figure(fignum)
        if clf:
            py.clf()
             
        py.scatter(pos_arr[:,0],pos_arr[:,1],
                        c=[i.type.color for i in self.cells],
                        s=50,
                        zorder=10,
                        *args,**kwargs)
        if ghosts:
            ghost_arr = self._get_ghost_pos_arr()
            py.scatter(ghost_arr[:,0],ghost_arr[:,1],
                        c = [i.original.type.color for i in self._ghosts],
                        s = 30,
                        zorder=10,
                        alpha = 0.3,
                        *args,**kwargs)
        py.axis('equal')

    def my_circle_scatter(self, axes, x_array, y_array, rad_array, col_array, **kwargs):
        for x, y, R, c in zip(x_array, y_array , rad_array, col_array):
            circle = py.Circle((x,y), radius=R, color = c, **kwargs)
            axes.add_patch(circle)
        return True

    def plot_sized_cells_old(self,clf=True,fignum=1,ghosts=False,*args, **kwargs):
        """ Plot the current configuration using circles"""
        self.logger.info('Plotting Sized Cells')
        pos_arr = self.get_pos_arr()
        rad_arr = self.get_radius_arr()
        col_arr = [i.type.color for i in self.cells]

        py.figure(fignum)
        if clf:
            py.clf()
           
        axes=py.axes()
        self.my_circle_scatter(axes,
                            pos_arr[:,0],
                            pos_arr[:,1],
                            rad_arr, col_arr, alpha=0.6,**kwargs)
                     
        if ghosts:
            ghost_arr = self._get_ghost_pos_arr()
            py.scatter(ghost_arr[:,0],ghost_arr[:,1],
                        c = [i.original.type.color for i in self._ghosts],
                        s = 30,
                        zorder=10,
                        alpha = 0.3,
                        *args,**kwargs)
        
        py.xlim((0,self.XSIZE))
        py.axis('equal')

    def plot_sized_cells(self,clf=True,fignum=1,ghosts=False,*args, **kwargs):
        """ Plot the current configuration using circles"""
        self.logger.info('Plotting Sized Cells')
        pos_arr = self.get_pos_arr()
        rad_arr = self.get_radius_arr()

        
        
        pos = self.get_pos_arr(force=True)
        pressure_arr = zeros_like(pos)
        

        #kdtree = self._get_kdtree(force=True)
        for i,j in self._get_npairs(): #kdtree.query_pairs(self.xi*1.0):
            
            force = self.force_func_celltypes(self.cells[i], self.cells[j] )
            pressure_arr[i] += fabs(force)
            pressure_arr[j] += fabs(force)
       
        pressure_arr = nan_to_num(pressure_arr)    

        #print "\n"
        #print pressure_arr
        #print "\n"
        
        cancer_cell_pressures = empty(len(self.cancer_cells))
        numero_cancer = 0
        numero_cell = 0
        for i in self.cells:
            if i.type.name == 'Cancer' :
                cancer_cell_pressures[numero_cancer]=norm(pressure_arr[numero_cell])/(3.141592*rad_arr[numero_cell]*rad_arr[numero_cell])
                numero_cancer =  numero_cancer + 1
            numero_cell =  numero_cell + 1

        #printing stress on file
        file=open(self.pressure_filename,'a')

        #factor is 4/3( E/(1-nu^2)) =  3/2 kPa 
        factor = 1.5

        for i in range(0,len(cancer_cell_pressures)):
            s = str(i) + ' ' + str(cancer_cell_pressures[i]*factor) +'\n'
            file.write(s)
            
        s = '\n'
        file.write(s)
        #s = str(numero_cancer) + ' ' + str(cancer_cell_pressures.mean()) +'\n'
            #file.write(s)
        #s = '\n'
        file.close
        

        
        if len(cancer_cell_pressures)>1 :
            cancer_cell_pressures = (cancer_cell_pressures-cancer_cell_pressures.min())/(cancer_cell_pressures.max()-cancer_cell_pressures.min())*0.9+0.1
            #print "\n"
            #print cancer_cell_pressures
            #print "\n"
        else :
            cancer_cell_pressures[0] = 0.5
       #print '\n'
        #print  cancer_cell_pressures
        #print '\n'
            


        col_arr = []
        numero_cancer = 0
        for i in self.cells:
            if i.type.name == 'Cancer' :
                rgb_color = cm.hot(1-cancer_cell_pressures[numero_cancer],1.0)
                col_arr.append(rgb_color)
                #print '\n'
                #print rgb_color , cancer_cell_forces[numero_cancer]
                #print '\n'
                numero_cancer =  numero_cancer + 1               
            else :
                col_arr.append(i.type.color)
        #print '\n'
        #print col_arr
        #print '\n'
        

        #file=open(self.screenshot_filename,'a')
        #for i in range(0, len(pos_arr)):
        #    s = self.cells[i].type.name + ' ' + str(pos_arr[i][0]) + ' ' + str(pos_arr[i][1]) + ' ' + str(rad_arr[i])  + ' ' + str(col_arr[i])  +'\n'
        #    file.write(s)
        #file.close
           

      


        py.figure(fignum)
        if clf:
            py.clf()
           
        axes=py.axes()
        self.my_circle_scatter(axes,
                            pos_arr[:,0],
                            pos_arr[:,1],
                            rad_arr, col_arr, alpha=0.6,**kwargs)
                     
        if ghosts:
            ghost_arr = self._get_ghost_pos_arr()
            py.scatter(ghost_arr[:,0],ghost_arr[:,1],
                        c = [i.original.type.color for i in self._ghosts],
                        s = 30,
                        zorder=10,
                        alpha = 0.3,
                        *args,**kwargs)
        
        py.xlim((0,self.XSIZE))
        py.axis('equal')



    def plot_links(self,clf=False,cutoff=None,fignum=1,ghosts=False,*args,**kwargs):
        """ Plot the links between cells """  
        self.logger.info('Plotting Links')
        if cutoff is None:
            cutoff = self.XSIZE/2.          
        py.figure(fignum)
        if clf:
            py.clf()
        
        #file=open(self.screenshot_filename,'a')


        for link in self.links:
            if link.C_10 > 0:
                #s = 'Link' + ' ' + str(link.one.pos[0]) + ' ' + str(link.one.pos[1]) + ' ' + str(link.two.pos[0])  +  ' ' + str(link.two.pos[1])  +'\n'
                #file.write(s)

                d12=link.one.pos-link.two.pos
                abs_d12=norm(d12)
                if abs_d12 < cutoff:
                    data = sp.array([ link.one.pos, link.two.pos ])
                    py.plot(data[:,0],data[:,1],
                                c=py.cm.jet( min(link.energy*30.,1.) ),
                                alpha=0.6,
                                *args, **kwargs )
        #file.close
    
    def _get_pairs(self):
        kdtree = self._get_kdtree(force=True)
        return kdtree.query_pairs(self.xi*1.0)

    def _get_cpairs(self,num=100):
        pos = self.get_pos_arr(force=True)
        ckdtree = self._get_ckdtree(force=False)
        ds,neighs = ckdtree.query(pos,num,distance_upper_bound=self.xi)

        pairs = set()
        N = len(neighs)
        for (i,j),k in sp.ndenumerate(neighs):
        #    if cmp(i,k) < 1:
        #        pairs.add((i,k))
        #    else:
        #        pairs.add((k,i))
            if k < N and (i,k) not in pairs and (k,i) not in pairs:
                pairs.add((i,k))

        return pairs

    def _get_npairs(self):
        if self.neighs is None:
            self.neighs = NeighborList([self.xi]*self.num_cells)
        
        self.neighs.update(self)

        return ((i,j) for i in range(self.num_cells) for j in self.neighs.get_neighbors(i) )


        

    @property
    def forces(self):
        """ get the forces between cells, as array, both from links
            and from the native force_func
        """
        self.logger.info('Computing forces')
        pos = self.get_pos_arr(force=True)

        force_arr = zeros_like(pos)

        for link in self.links:
            force = link.force
            force_arr[link.one.index] += force
            force_arr[link.two.index] -= force


        #kdtree = self._get_kdtree(force=True)
        for i,j in self._get_npairs(): #kdtree.query_pairs(self.xi*1.0):
            
            force = self.force_func_celltypes(self.cells[i], self.cells[j] )
            #disp = self.cells[i].pos - self.cells[j].pos
            #L = norm(disp)
            #force = 2 * self.a**4 * ( 2 * self.xi**2 - 3 * self.xi * L + L**2 )/( self.xi**2 * L**6 ) * disp
            force_arr[i] += force
            force_arr[j] -= force

        return nan_to_num(force_arr)


    def force_func(self,cell1,cell2):
        """ the native force function between two positions """
        x1 = cell1.pos
        x2 = cell2.pos
        disp = x1 - x2
        mod_disp = norm(disp)
        force = 2 * self.a**4 * ( 2 * self.xi**2 - 3 * self.xi * mod_disp + mod_disp**2 )/( self.xi**2 * mod_disp**6 ) * disp
        
        return force

    def force_func2(self,cell1,cell2):
        """ the native force function between two positions, second attempt """
        x1 = cell1.pos
        x2 = cell2.pos
        r1 = cell1.radius
        r2 = cell2.radius
        disp = x1 - x2
        mod_disp = norm(disp)
        a1=self.a*(r1+r2)
        xi1=self.xi*(r1+r2)
        force = 2 * a1**4 * ( 2 * xi1**2 - 3 * xi1 * mod_disp + mod_disp**2 )/( xi1**2 * mod_disp**6 ) * disp
        
        return force

    def force_func_hertz(self,cell1,cell2):
        """ the Hertz force between two cells """
        x1 = cell1.pos
        x2 = cell2.pos
        r1 = cell1.radius
        r2 = cell2.radius
        disp = x1 - x2
        mod_disp = norm(disp)
        delta=(r1+r2)-mod_disp
        if delta > 0.0:
            force = self.a*delta**1.5*disp/mod_disp
        else:
            force= 0.0

        return force

    def force_func_celltypes_old(self,cell1,cell2):
        """ Try to case out the cell types """

        x1 = cell1.pos
        x2 = cell2.pos

        #use the Cython dispfunc
        disp = disp_func(x1,x2,self.XSIZE)
        mod_disp = norm(disp)
        force = 0.0

        if cell1.type==self.basal and cell2.type==self.basal:
            #We have two basal cells
            force = 0.0
        elif cell1.type==self.basal or cell2.type==self.basal:
            #We have one basal cell
            if mod_disp <= self.basalcutoff:
                oldexpr = '2 * self.basalstrength**4 * ( 2 * self.basalcutoff**2 - 3 * self.basalcutoff * mod_disp + mod_disp**2 )/( self.basalcutoff**2 * mod_disp**6 ) * disp'
                
                basalstrength = self.basalstrength
                basalcutoff = self.basalcutoff
                forcestr = '2 * basalstrength**4 * ( 2 * basalcutoff**2 - 3 * basalcutoff * mod_disp + mod_disp**2 )/( basalcutoff**2 * mod_disp**6 ) * disp'
                force = ne.evaluate(forcestr)
        else:
            #We have some other situation

            r1 = cell1.radius
            r2 = cell2.radius
            delta=(r1+r2)-mod_disp

            if delta > 0:
                a = self.a
                oldexp = 'sqrt(r1*r2/(r1+r2)) * self.a * delta**1.5*disp/mod_disp'
                forcestr = 'sqrt(r1*r2/(r1+r2)) * a * delta**1.5*disp/mod_disp'
                force = ne.evaluate(forcestr)
                #print 'force', force

        return force

    def force_func_celltypes(self,cell1,cell2):
        """ Try to case out the cell types """

        x1 = cell1.pos
        x2 = cell2.pos

        #use the Cython dispfunc
        disp = disp_func(x1,x2,self.XSIZE)
        mod_disp = norm(disp)
        force = 0.0

        if cell1.type==self.basal and cell2.type==self.basal:
            #We have two basal cells
            force = 0.0
        #elif cell1.type==self.basal or cell2.type==self.basal:
            #We have one basal cell
         #   if mod_disp <= self.basalcutoff:
         #       oldexpr = '2 * self.basalstrength**4 * ( 2 * self.basalcutoff**2 - 3 * self.basalcutoff * mod_disp + mod_disp**2 )/( self.basalcutoff**2 * mod_disp**6 ) * disp'
                
          #      basalstrength = self.basalstrength
          #      basalcutoff = self.basalcutoff
          #      forcestr = '2 * basalstrength**4 * ( 2 * basalcutoff**2 - 3 * basalcutoff * mod_disp + mod_disp**2 )/( basalcutoff**2 * mod_disp**6 ) * disp'
          #      force = ne.evaluate(forcestr)
        else:
            #We have some other situation

            r1 = cell1.radius
            r2 = cell2.radius
            
            min_radius = min(r1,r2)

            renormalized_r = r1*r2/(r1+r2)

            delta=(r1+r2)-mod_disp
            
            if delta > 0:
                omega = pow(delta/renormalized_r,1.5)
                a = self.a
                forcestr = 'sqrt(renormalized_r) * a * delta**1.5*(1 + 1.15*omega**0.34 +9.5*omega + 9.288*omega**2)/(1+2.3*omega)*disp/mod_disp'
                force = ne.evaluate(forcestr)
                
        return force


    def force_func_celltypes_cython(self,cell1,cell2):
        """ Try to case out the cell types """

        x1 = cell1.pos
        x2 = cell2.pos


        if cell1.type==self.basal and cell2.type==self.basal:
            #We have two basal cells
            force = 0.0
        elif cell1.type==self.basal or cell2.type==self.basal:
            #We have one basal cell
            force = force_func_basal(x1,x2,self.basalstrength,self.XSIZE)
        else:
            #We have some other situation
            r1 = cell1.radius
            r2 = cell2.radius
            force = force_func_hertz(x1,x2,r1,r2,self.a,self.XSIZE)

        return force


    @property
    def energy(self):
        """ get the energy of the current configuration """
        tot_energy = 0
        for link in self.links:
            tot_energy += link.energy
        return tot_energy

    def fire(self):
        """ Do a fire relaxation """

        #load params
        fmax = self.config['fmax']
        Nmin = self.config['Nmin']
        finc = self.config['finc']
        fdec = self.config['fdec']
        alphastart = self.config['alphastart']
        fa = self.config['fa']
        deltatmax = self.config['deltatmax']
        maxsteps = self.config['maxsteps']


        alpha = alphastart
        deltat = 0.1

        pos = self.get_pos_arr(force=True)
        v = sp.zeros_like(pos)
        self._update_vel(v)

        v = self._get_vel_arr()

        steps_since_negative = 0

        def norm_arr_old(vec):
            return sp.sqrt(sp.sum(vec**2,1))
        def unitize_arr_old(vec):
            return nan_to_num(((vec.T)/norm_arr(vec)).T)

        norm_arr = norm

        forces = nan_to_num(sp.array([ [sp.inf,sp.inf]]))

        step_num = 0

        self.logger.info("Beginning FIRE Relaxation -- fmax={}".format(fmax))
        
        maxdpos = 100000.0
                 

        while max(norm_arr(forces)) > fmax and step_num < maxsteps:
            forces = self.forces

            self.logger.debug("Computed forces: {forces}".format(forces=pprint.pformat(forces)))

            power = sp.vdot(forces,v)
            self.logger.info("Step: {}, max_force: {}, power: {}".format(step_num,
                                        max(norm_arr(forces)),
                                         power))
            
            #DEBUG PRINTING
            #print "Step: {}, max_force: {}, power: {}, deltat: {}".format(step_num,
            #                            max(norm_arr(forces)),
            #                             power, deltat)

            v = nan_to_num( (1.0 - alpha)*v + alpha*(norm_arr(v)*unitize_arr(forces).T).T )

            if power>0.:
                if steps_since_negative > Nmin:
                    deltat = min(deltat * finc, deltatmax)
                    alpha = alpha*fa
                steps_since_negative += 1

            else:
                steps_since_negative = 0

                deltat = deltat * fdec
                v *= 0.
                alpha = alphastart

            v += forces*deltat
            pos += v*deltat

            self._update_pos(pos)
            step_num += 1

            #maxdpos = max(norm_arr(v*deltat))
            #DEBUG PRINTING
            #print "Maximum position change = {}".format(maxdpos)
            #DEBUG_PLOT
            #self.plot_sized_cells()
            #self.plot_links()
            #self.plot_forces()
            #py.draw()

        self._update_pos(pos)
        self._update_vel(v)
        self.logger.info("Relaxation finished...")
    

    def save(self,filename):
        self.logger.info("SAVING state to {}".format(filename))
        with open(filename,'w') as f:
            pickle.dump( (self.config, self.cells, self.links, self._ghosts, self.T ), f )



    def vmd_out(self,filename):
        """ Write a VMD compatible file to filename """
        with open(filename,'w') as f:
            
            positions = self.get_pos_arr(force=True)

            formatstring  = "{color} {x} {y} {z}\n"

            for ind,row in enumerate(positions):
                f.write(formatstring.format(x=row[0], y=row[1], z=0, color=self.cells[ind].type.type_ind))


    def plot_forces(self,factor=5):
        X,Y = self.get_pos_arr().T
        FX,FY = self.forces.T

        py.quiver(X,Y,FX,FY,scale=factor)


    #Some code for ASE neighborlist functionality
    def get_positions(self):
        return sp.hstack(( self.get_pos_arr(), sp.zeros((self.num_cells,1)) ) )

    def get_pbc(self):
        return sp.array([True,False,False])
    
    def get_cell(self):
        return sp.array([[self.XSIZE,0,0],[0,self.YSIZE,0],[0,0,1]])

    def __len__(self):
        return self.num_cells



def load_from_file(filename):
    with open(filename,'r') as f:
        config, cells, links, ghosts, T = pickle.load(f)
    Q = CancerSim(config)
    Q.cells = cells
    Q.ghosts = ghosts
    Q.T = T
    Q.links = links
    Q.cancer_cells = [cell for cell in cells if cell.type.name == "Cancer"]
    Q.num_cells = len(Q.cells)
    return Q


if __name__ == "__main__":
    Q = CancerSim()
    Q._triang_lattice()
    Q.delaunay()
    Q._freeze_links()

    Q.add_cancer_cell([XSIZE/2.,YSIZE/2 + 3],1)

    Q.plot_cells()


    self = Q


"""
TODO:  have links know about periodic boundary conditions (maybe)
        freeze links (DONE)
        Ghost cells need update method.  (DONE)
        fire relaxation (DONE)
        set and divide cancer cells (DONE)
        long range forces (DONE)
            
        cache the link calcs
        cache the KDTree calcs?
        allow more transparent custimization
        expose CellTypes
        use logging module
"""
