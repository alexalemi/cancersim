import scipy as sp

from math import sqrt
from scipy import dot

from scipy.linalg import norm as builtinnorm

########################################################
### Helper functions ###################################
########################################################

import numexpr as ne
from numpy import newaxis
from scipy import nan_to_num

def disp_func(x1,x2,XSIZE):
    disp1 = ne.evaluate('x1 - x2')
    norm1 = norm(disp1)

    if norm1 < 3:
    	return disp1

    disp2 = x1 + XSIZE - x2
    disp3 = x1 - XSIZE - x2

    norm1 = norm(disp1)
    norm2 = norm(disp2)
    norm3 = norm(disp3)

    if norm1 <= norm2 and norm1 <= norm3:
        return disp1
    elif norm2 <= norm1 and norm2 <= norm3:
        return disp2
    else:
        return disp3


def normold(vec):
    """ Normalize a vector """
    return sp.sqrt(sp.sum(vec**2))


def norm(vec):
    temp = ne.evaluate('sum(vec**2)')
    return ne.evaluate('sqrt(temp)')


def normnp(vec):
	""" Second norm """
	return sqrt(dot(vec,vec.conj()))

def norm3(vec):
	""" Third norm """
	return builtinnorm(vec)

def unitize(vec):
    """ Construct a unit vector from a vector """
    return vec/norm(vec)

def unitize_arr(vec):
    norms = norm(vec)[:,newaxis]
    return nan_to_num(ne.evaluate('vec/norms'))
