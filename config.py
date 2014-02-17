from __future__ import division

""" 
Configuration
this sets the parameters for the simulation


As for units, we are working in units where our natural length scale is 1e-5 m or 10 um

and so that 
A = 4/3( E/(1-nu^2)) = 4/3 ( (1 kPa) / ( 1 - (1/3)^2 ) = 3/2 kPa = 1 

And with a natural adhesion between cells of 200 uN/m, we have the factor below
"""









config = {
	#General params
        'XSIZE': 60, 
	'YSIZE': 60, 
	'basal_height': 37.5,
        'basal_wavenumber': 2.0, 
        'basal_amplitude': 7.5, 
	'seed': 1307,
	'force_magnitude': 1.0,  

        #Pressure filename
        'pressure_filename' : 'pressure.dat',
        'cancer_evolution_filename' : 'cancer_evolution.dat',

	'force_magnitude_basal': 0.10,
	'force_cutoff': 2.0,
	'force_cutoff_basal': 0.5,
	'jiggle_sigma': 0.1,

	#First Cancer Cell
	'first_cancer_cell_yoffset': 1.,
        'first_cancer_cell_xoffset': 0.,
	'first_cancer_cell_radius': 5.,

	#CELL PARAMETERS
	'cancer_cell_params': { 
		'name':'Cancer',
		'type_ind':0,
		'C_10': 0.0,
                'C_11': 0.0,
		'L': 2.,
		'color': 'y',
		'maxstretch': 1.2,
	},

	'epidermal_cell_params': {
		'name': 'Epidermal',
		'type_ind': 1,
		'C_10': 0.0075,
                'C_11': 0.2,
		'L': 2.0,
		'color': 'b',
		'maxstretch': 1.2,
	},

	'basal_cell_params': {
		'name': 'Basal',
		'type_ind': 2,
		'C_10': 0.1,
                'C_11': 0.0,
		'L': 1.0,
		'color': 'k',
		'maxstretch': 1.2,
	},

	'dermal_cell_params': {
		'name': 'tDermal',
		'type_ind': 3,
		'C_10': 0.0075,
                'C_11': 0.2,
		'L': 6.0,
		'color': 'g',
		'maxstretch': 1.2,
	},

        'stratum_corneum_cell_params': {
		'name': 'Corneum',
		'type_ind': 4,
		'C_10': 1.0,
                'C_11': 0.0,
		'L': 2.0,
		'color': 'r',
		'maxstretch': 1.3,
	},


	#Parameters for FIRE
	'fmax': 0.05,
	'Nmin': 5.,
	'finc': 1.1,
	'fdec': 0.5,
	'alphastart': 0.1,
	'fa': 0.99,
	'deltatmax': 10.,
        'maxsteps': 2*10**3,
}
