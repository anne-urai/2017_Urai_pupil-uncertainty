#!/usr/bin/env python

__doc__ = """GLM fitting with weighted data points

Copyright (C) 2014 Ingo Fruend

This code reproduces the analyses in the paper

    Fruend, Wichmann, Macke (2014): Quantifying the effect of inter-trial dependence on perceptual decisions. J Vis, 14(7): 9.


    Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

    If you use the Software for your own research, cite the paper.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

import numpy as np
import sys

def optimize_w ( X, y, q, w, niter=5, stop=1e-5, lm=0.1, verbose=False ):
    """Optimize the w parameters of the model with weighted data

    :Parameters:
        *X*     design matrix
        *y*     responses
        *q*     state probabilities (weights)
        *w*     starting weights

    :Optional parameters:
        *niter* number of iterations
        *stop*  stopping criterion
        *lm*    regularization
        *verbose* obvious
    """
    w,l_,conv = weighted_glm ( X, y, q, w, niter=niter, stop=stop, lm=lm )
    if verbose:
        print "Converged" if conv else "Not converged"
    return np.ravel(np.array(w))

def logistic ( x ):
    """Logistic function"""
    return 1./(1+np.exp(-x))

def weighted_glm ( X, y, q, w, niter=5, stop=1e-5, lm=0.1 ):
    """The actual optimization

    Parameters: see optimize_w
    """
    X = np.matrix ( X )
    w = np.matrix ( w.reshape((-1,1)) )
    y = np.matrix ( y.reshape((-1,1)) )
    q = np.matrix ( q.reshape((-1,1)) )
    eta = logistic ( X*w )
    l = np.sum(q.A*(y.A*np.log(eta.A)+(1-y.A)*np.log(1-eta.A))) - w.T*w

    for iteration in xrange ( niter ):
        z = q.A*(y-eta).A
        grad_l = X.T*z - 2*lm*w
        w_ii = -q.A*eta.A*(1-eta.A)
        H = (X.A*w_ii.reshape ( (-1,1) )).T * X
        H -= 2*lm*np.eye ( len(w) )

        dw = np.linalg.lstsq ( H, grad_l )[0]
        e = np.sqrt ( np.sum ( (dw).A**2 ) )

        w -= 0.9*dw
        eta = logistic ( X*w )

        if e < stop:
            break

    l_ = np.sum(q.A*(y.A*np.log(eta.A)+(1-y.A)*np.log(1-eta.A)))
    wnorm = - w.T*w

    if not l_-wnorm>=l:
        conv = 0.
    else:
        conv = 1.

    return w,l_,conv

