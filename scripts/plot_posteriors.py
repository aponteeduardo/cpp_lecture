#! /usr/bin/env python

# aponteeduardo@gmail.com
# copyright (C) 2015


''''

Plot the posterior distributions.

'''

import os

import numpy as np
import scipy.io as scio
import pylab as pyl


def posteriors_fname():
    ''' Return name of the file with the posteriors. '''

    fname = 'posteriors.mat'
    dname = '../dump/posteriors'

    fname = os.path.join(dname, fname)

    return fname


if __name__ == '__main__':
    
    fname = posteriors_fname()
    post = scio.loadmat(fname)
    post = post['post']
    
    theta = post['theta'][0, 0]

    alpha = np.arctan(theta[0, :])/np.pi + 0.5
    pyl.plot(alpha)

    print(np.mean(alpha))

    pyl.figure()

    ax = pyl.hist(alpha, color='r', bins=50)
#    pyl.title(r'$\text{Distribution} \alpha$', size=30)
     
    beta = np.exp(theta[1, :])
    print np.mean(beta)

    pyl.figure()
    pyl.plot(beta)

    pyl.figure()
    
    pyl.hist(beta, color='g', bins=50)


    pyl.show()
    pass

