#!/usr/bin/env python

'''
hps.py

Calculates the number of times an air particle would collide with a pollen particle every
second. This only counts first interactions, and would only count an interaction that led
to a double bounce once.
'''

import numpy as np
import numpy.linalg as la
from scipy.spatial import ConvexHull
import itertools
import sys

def sphere2cart(s, stretch = 10):
   return stretch * np.array([np.sin(s[0]) * np.cos(s[1]), np.sin(s[0]) * np.sin(s[1]), np.cos(s[0])])

def hps(MESH_FN = '../meshes/in.obj', NUM_SAMP = 1000):
    '''
    Calculate hits per second

    MESH_FN:    filename of the obj file to calculate mass diffusivity for

    NUM_SAMP:   number of theta, phi to integrate over; note that actual 
                number of computations is NUM_SAMP ^ 2
    '''
    points = []
    with open(MESH_FN, 'r') as fi:
            print 'Start processing obj file'
            lines = fi.readlines()

            for i in xrange(len(lines)):
                    if lines[i][0] == 'v':
                            points.append([float(rep) for rep in lines[i][2:-1].split(' ')])
                    else:	# obj files do verticies then other stuff
                            break

            M = np.asmatrix(points).T	# got all the verticies
            print 'End processing obj file'

            # Generate all the normal vectors we are going to integrate over
            print 'Start generating all normal vectors'
            phis = np.linspace(0, 2 * np.pi, NUM_SAMP)
            thetas = np.linspace(0, 2 * np.pi, NUM_SAMP)

            cpi = itertools.product(*[phis, thetas])	# combination product iterator
            ncpi = map(np.array, cpi)
            cncpi = map(sphere2cart, ncpi)
            N = np.asmatrix(cncpi).T

            print 'Shape N: ' + str(np.shape(N))

            # Get basis of the normal vector
            BT_1 = np.matrix('0 0 1; 0 0 0; -1 0 0') # basis transform
            BT_2 = np.matrix('0 0 0; 0 0 1; 0 -1 0')

            # assume N is normal vec
            print 'Start getting all orthonormal basis'
            bas_1 = BT_1 * N
            bas_2 = BT_2 * N

            print 'Shape bas_1: ' + str(np.shape(bas_1))

            col_it = (np.column_stack([bas_1[:,i], bas_2[:,i]]) for i in xrange(np.shape(N)[1]))
            obas = []

            for mat in col_it:
                    obas.append(la.qr(mat)[0])

            print 'End getting all orthonormal basis'

            # Assume M is a set of coordinates
            print 'Start getting all convex hulls'
            rep = 0.0
            projs_it = (bas_mat.T * M for bas_mat in obas)


            c = sum((ConvexHull(proj.T).area for proj in projs_it))
            print 'End getting all convex hulls'

            return c

            print 'Final answer: ' + str(c)

if __name__ == '__main__':
    if len(sys.argv) == 1:
        print 'Final answer: ' + str(hps())
    elif len(sys.argv) == 2:
        if sys.argv[1].isdigit():
            print 'Final answer: ' + str(hps(NUM_SAMP=int(sys.argv[1])))
        else:
            print 'Final answer: ' + str(hps(MESH_FN = sys.argv[1]))
    elif len(sys.argv == 3):
        print 'Final answer: ' + str(hps(MESH_FN = sys.argv[1], NUM_SAMP=int(sys.argv[2])))
    else:
        print 'Too many args'
