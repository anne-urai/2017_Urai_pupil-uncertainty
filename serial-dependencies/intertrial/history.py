#!/usr/bin/env python

import numpy as np
import copy,sys

__doc__ = """This module defines the basic functionality to define a dataset with history features in it.
The central component of this module is the DataSet class. It is meant to be a base class for data sets
that can be used in the analysis of history features.

One of the central components of the DataSet class is its ability to generate random datasets -- permutations
of the original dataset or random resamples.


Copyright (C) 2014 Ingo Fruend

This code reproduces the analyses in the paper

    Fruend, Wichmann, Macke (2014): Quantifying the effect of inter-trial dependence on perceptual decisions. J Vis, 14(7): 9.


    Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

    If you use the Software for your own research, cite the paper.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

__all__ = ["DataSet","history_impulses","history_features","get_w0"]

def modifiable_property(function):
    keys = 'fget', 'fset', 'fdel'
    func_locals = {'doc':function.__doc__}
    def probe_func(frame, event, arg):
        if event == 'return':
            locals = frame.f_locals
            func_locals.update(dict((k, locals.get(k)) for k in keys))
            sys.settrace(None)
        return probe_func
    sys.settrace(probe_func)
    function()
    return property(**func_locals)

#################### A data set structure ####################

class DataSet ( object ):
    """Base class for all datasets"""
    def __init__ ( self, impulse_responses ):
        """
        :Parameters:
            *impulse_responses*
                A matrix of shape (nlags,nhist) with one history feature per column. If
                this is a vector, it will be interpreted as the impulse response for a
                single history feature

        :Example:
        >>> d = DataSet ( [[1.,1.],[0.,1.]] )
        >>> d = DataSet ( [1.,0.] )
        >>> d = DataSet ( [[[1.,1.],[0.,1.]]] )
        Traceback (most recent call last):
        ...
        ValueError: Don't know how to handle an impulse response matrix of this dimension
        """
        if impulse_responses is None:
            impulse_responses = []
        impulse_responses = np.array ( impulse_responses )
        if len ( impulse_responses.shape )==1:
            impulse_responses = np.reshape ( impulse_responses, (-1,1) )
        elif len ( impulse_responses.shape ) > 2:
            raise ValueError("Don't know how to handle an impulse response matrix of this dimension")
        self.h = gram_schmidt ( impulse_responses )

        self.__fname = "no file"

        self.rng = np.random.mtrand.RandomState ()

    # Some "virtual" methods that need to be implemented by the user
    def bootstrap(self):
        """return a bootstrap dataset

        Derived classes might not need to overwrite this one
        """
        i = self.rng.randint ( len(self.r), size=len(self.r) )
        return self.r[i],self.X[i,:]
    def permutation ( self ):
        """Generate response vector and design matrix  of a random permutation of the data set"""
        raise NotImplementedError
    @property
    def X ( self ):
        """Return the design matrix"""
        raise NotImplementedError
    @property
    def r ( self ):
        """Return the response vector"""
        raise NotImplementedError

    def getsummary ( self, condition=0 ):
        """Return a datasummary for one condition"""
        raise NotImplementedError

    @property
    def hf0 ( self ):
        return 2

    @property
    def th_features ( self ):
        """Features to be thresholded"""
        return []

    # These do not need to be implemented by the user
    def gethistorykernel ( self, w, alpha=1. ):
        """Return the history kernel associated with the respective weights

        :Parameters:
            *w*
                weights for the different history basis functions
            *alpha*
                model coefficient associated with the stimulus

        :Example:
        >>> impulse_responses = [[1.,1.,1.],[0.,1.,1.],[0.,0.,1.],[0.,1.,0.]]
        >>> d = DataSet ( impulse_responses )
        >>> d.gethistorykernel ( [1.,0.,0.] )
        array([ 1.,  0.,  0.,  0.])
        """
        w = np.matrix ( w ).T
        return np.ravel ( np.array ( np.matrix ( self.h )*w ) )/alpha

    @property
    def hlen ( self ):
        return self.h.shape[1]

    @modifiable_property
    def fname ():
        def fget ( self ):
            return self.__fname
        def fset ( self, v ):
            self.__fname=v

#################### Gram Schmidt ####################

def proj ( v, u ):
    """projection of v on u

    :Example:
    >>> v = [1.,2.]
    >>> u = [0.,1.]
    >>> proj ( v, u )
    array([ 0.,  2.])
    """
    w = np.dot ( v, u )/np.dot ( u,u )
    return w*np.array(u)

def vector_norm ( v ):
    """return the vector norm

    :Example:
    >>> v = [3.,4.]
    >>> vector_norm ( v )
    5.0
    """
    return np.sqrt ( np.dot ( v, v ) )

def gram_schmidt ( X ):
    """Perform Gram-Schmidt orthogonalization of the columns

    :Example:
    >>> X = [[1.,1.,1.],[0.,1.,1.],[0.,0.,1.],[1.,0.,0.]]
    >>> E = gram_schmidt ( X )
    >>> np.dot ( E.T, E )
    array([[ 1.,  0.,  0.],
           [ 0.,  1.,  0.],
           [ 0.,  0.,  1.]])
    """
    X = np.array ( X )
    assert len ( X.shape ) == 2

    U = [X[:,0]]
    E = [X[:,0]/vector_norm ( X[:,0] )]

    for i in range ( 1, X.shape[1] ):
        u_new = X[:,i].copy()
        v = X[:,i]
        for u in U:
            u_new -= proj ( v, u )
        U.append ( u_new )
        E.append ( u_new/vector_norm ( u_new ) )

    return np.transpose(E)

#################### Dealing with history features ####################

# sys.stderr.write ( "WARNING USING ONLY 1Back feature!\n" )
# def history_impulses ( eta=[-1], nlags=7 ):

# original
# def history_impulses ( eta=[-1,0.5,0.25], nlags=7 ):

# AEU: change this to independent lags
# def history_impulses ( eta=[-1,-2,-3,-4,-5,-6,-7], nlags=7 ):
def history_impulses ( eta=[-1,0.5,0.25], nlags=7 ):

    """Determine impulse responses of history features

    :Parameters:
        *eta*
            a sequence of scaling factors or lags for the different history features
            Negative values will be considered as lags, positive values as scaling
            factors

    :Example:
    >>> history_impulses ( [-1,-2,0.5], 4 )
    array([[ 1.    ,  0.    ,  0.5   ],
           [ 0.    ,  1.    ,  0.25  ],
           [ 0.    ,  0.    ,  0.125 ],
           [ 0.    ,  0.    ,  0.0625]])
    >>> history_impulses ( [-1.1] )
    Traceback (most recent call last):
    ...
    AssertionError: negative history parameters will be interpreted as lags and should be integers
    """
    impulse_responses = np.zeros ( (nlags,len(eta) ) )
    for i,e in enumerate ( eta ):
        if e>0:
            impulse_responses[:,i] = e**np.arange ( 1, nlags+1 )
        elif e<0:
            assert (e%1) == 0, "negative history parameters will be interpreted as lags and should be integers"
            impulse_responses[-1-e,i] = 1.
    return impulse_responses

def history_features ( h, r_or_z ):
    """Determine the values of a history feature

    :Parameters:
        *h*
            The impulse responses for all history features. This should be a matrix containing the impulse
            responses in the colums
        *r_or_z*
            a vector of the data that should go into the history (i.e. responses r or stimulus identities z).

    :Returns:
        a matrix with one column for every column in "imp" and one row for every entry in "r_or_z"

    :Example:
    >>> z = [1,-1,1,1,-1,-1,-1,-1,1]
    >>> h = history_impulses ( [-1,0.5] )
    >>> f = history_features ( h, z )
    >>> np.c_[z,f]
    array([[ 1.       ,  0.       ,  0.       ],
           [-1.       ,  1.       ,  0.5      ],
           [ 1.       , -1.       , -0.25     ],
           [ 1.       ,  1.       ,  0.375    ],
           [-1.       ,  1.       ,  0.6875   ],
           [-1.       , -1.       , -0.15625  ],
           [-1.       , -1.       , -0.578125 ],
           [-1.       , -1.       , -0.7890625],
           [ 1.       , -1.       , -0.8984375]])
    >>> h = np.array([[]])
    >>> history_features ( h, z )
    """
    if 0 in h.shape:
        return None

    assert len ( h.shape ) == 2, "h should be 2d"
    r_or_z = np.array ( r_or_z, 'd' )
    hf = np.zeros ( (len(r_or_z), h.shape[1]), 'd' )
    n,m = hf.shape

    for i in xrange ( m ):
        hf[1:,i] = np.convolve ( h[:,i], r_or_z, 'full' )[0:n-1]

    # added: orthogonalize design matrix
    # hf = gram_schmidt (hf)

    return hf

def history_features_stim (h, r_or_z, d):
    """Determine the values of a history feature

    :Parameters:
        *d* Difficulty of each trial, if 0 there will be no stimulus weight assigned
        *h*
            The impulse responses for all history features. This should be a matrix containing the impulse
            responses in the colums
        *r_or_z*
            a vector of the data that should go into the history (i.e. responses r or stimulus identities z).

    :Returns:
        a matrix with one column for every column in "imp" and one row for every entry in "r_or_z"
    """

    if 0 in h.shape:
        return None

    assert len ( h.shape ) == 2, "h should be 2d"
    r_or_z = np.array ( r_or_z, 'd' )
    hf = np.zeros ( (len(r_or_z), h.shape[1]), 'd' )
    n,m = hf.shape

    for i in xrange ( m ):
        # see slack message Anke 15 July 2017: when there is no stimulus evidence, also set this history feature to zero
        for j in xrange ( len(r_or_z)):
            if d[j] == 0:
                r_or_z[j] = 0
        hf[1:,i] = np.convolve ( h[:,i], r_or_z, 'full' )[0:n-1]

    # added: orthogonalize design matrix
    # hf = gram_schmidt (hf)

    return hf

def get_code ( x, code_out, code_in ):
    """Convert from one coding to another

    :Parameters:
        *x*
            input to be coded
        *code_out*
            desired code of the output
        *code_in*
            code of the input

    :Example:
    >>> code_out = [-1,1]
    >>> code_in  = [0,1]
    >>> get_code ( 1., code_out, code_in )
    1
    >>> get_code ( 0., code_out, code_in )
    -1
    >>> get_code ( 0., code_in, code_in )
    0
    >>> get_code ( 0., code_in, code_out )
    Traceback (most recent call last):
    ...
    ValueError: list.index(x): x not in list

    """
    if isinstance ( code_in, (list,tuple) ):
        return code_out[code_in.index(x)]
    else:
        return code_out[code_in.tolist ().index(x)]

################## starting values for optimization ##################

def get_w0 ( d ):
    """Determine starting values for stimulus related parameters

    :Parameters:
        *d*
            a matrix of the same form as returned by DataSet.getsummary(). The
            first column is the (signed) stimulus levels, the second column is
            the fraction of responses that are to be modeled, and the third
            column is the number of trials at that particular stimulus level"""
    x = d[:,0]
    k = d[:,1]
    n = d[:,2]

    assert not (n==0).any()

    # Avoid zero and one by adding half a trial in the respective cases
    k[k==0] += 0.5
    n[k>=n] += 0.5

    p = k/n
    mu = np.log ( p/(1-p) )
    X = np.ones ( (len(x),2), 'd' )
    X[:,1] = x

    # Now use a Taylor expansion of the model around the inflection point x0 = -w0/w1:
    # f(x) \approx 1/2 + 1/4 w1 ( x + w0/w1 )
    a,b = np.linalg.lstsq ( X, mu )[0]
    w0,w1 = 4*(a-0.5),4*b

    return np.array( [w0,w1] )

if __name__ == "__main__":
    import doctest
    doctest.testmod ()
