#!/usr/bin/env python

import numpy as np
import scipy.optimize as opt
import model,threshold
import sys

__doc__ = """Functions to perform the monte carlo simulations for determining statistical significance

The central function of this module is mcsimulation, that takes two arguments that might need some explanation

generator will typically be one of the DataSet methods permutation() and bootstrap()
In general it should return a response vector and a design matrix

collector, should take a fitted model and return the values that should be stored from that model. Three examples
are in this module. The KernelCollector returns the history kernels of the fitted model -- this will typically be
needed for bootstrap, the LikelihoodCollector returns only the log-likelihood and the Likelihood_and_Slope_Collector
returns log-likelihood and slope. The latter two collectors will typically be used for the permutation tests.


Copyright (C) 2014 Ingo Fruend

This code reproduces the analyses in the paper

    Fruend, Wichmann, Macke (2014): Quantifying the effect of inter-trial dependence on perceptual decisions. J Vis, 14(7): 9.


    Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

    If you use the Software for your own research, cite the paper.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

#################### Some potential collector objects ####################

class KernelCollector ( object ):
    def __init__ ( self, impulse_responses, starting_index, slopeindex=1 ):
        """A collector object for the history_kernels

        :Parameters:
            *impulse_responses*
                an array containing the (orthogonal) impulse responses in the columns
                This will usually be provided by a DataSet d as d.h
            *starting_index*
                index of the first history feature

        :Returns:
            all kernels concatenated into one long array

        :Example:
        >>> kcollect = KernelCollector ( np.eye(3), 3 )
        >>> class fake_model: pass
        >>> M = fake_model()
        >>> M.w = np.arange ( 9 )
        >>> kcollect ( M )
        array([ 3.,  4.,  5.,  6.,  7.,  8.])
        """
        assert impulse_responses.ndim == 2
        self.h  = impulse_responses
        self.i0 = starting_index
        self.slopeindex = slopeindex

    def __call__ ( self, fitted_model ):
        out = []
        w = fitted_model.w[self.i0:]
        nbase = self.h.shape[1]
        assert len(w)%nbase == 0

        if getattr ( self.slopeindex, "__iter__", False ):
            alpha = fitted_model.w[self.slopeindex].mean()
        else:
            alpha = fitted_model.w[self.slopeindex]

        for i in xrange ( 0, len(w), nbase ):
            out.append ( np.dot ( self.h, w[i:i+nbase] )/alpha )

        if len(out)==0:
            return np.array([])
        else:
            return np.concatenate ( out )

class Kernel_and_Slope_Collector ( KernelCollector ):
    def __init__ ( self, impulse_responses, starting_index, slopeindex=1 ):
        """A collector object for both, history kernels and slope

        :Parameters:
            *impulse_responses*
                an array containing the (orthogonal) impulse responses in the columns
                This will usually be provided by a DataSet d as d.h
            *starting_index*
                index of the first history feature
            *slopeindex*
                index of the slope parameter

        :Returns:
            all kernels concatenated into one long array
        """
        KernelCollector.__init__ ( self, impulse_responses, starting_index, slopeindex )
    def __call__ ( self, fitted_model ):
        out = KernelCollector.__call__ ( self, fitted_model )

        if getattr ( self.slopeindex, "__iter__", False ):
            alpha = fitted_model.w[self.slopeindex].mean()
        else:
            alpha = fitted_model.w[self.slopeindex]

        lm = fitted_model.pi[0]

        return np.concatenate ( (out,[alpha,lm]) )

def LikelihoodCollector (fitted_model):
    """only return the likelihood of the fitted model"""
    return fitted_model.loglikelihood

def Likelihood_and_Slope_Collector (fitted_model):
    """Collect likelihood and slope

    The output will be (loglikelihood,w(stimulus),slope_of_the_full_model_at_inflection)
    
    This assumes that the first two columns of the design matrix refer to the constant term and a term
    corresponding to the stimulus level

    :Example:
    >>> class fake_model: pass
    >>> M = fake_model ()
    >>> M.w  = np.ones( 2,'d' )
    >>> M.pi = np.array([.1,.1,.8])
    >>> M.loglikelihood = -920.
    >>> Likelihood_and_Slope_Collector ( M )
    array([ -9.20000000e+02,   1.00000000e+00,   2.00000000e-01])
    """
    out = np.zeros ( 3, 'd' )
    out[0] = fitted_model.loglikelihood
    out[1] = fitted_model.w[1]
    out[2] = fitted_model.pi[2]*0.25*fitted_model.w[1]

    return out

class EvaluationCollector ( object ):
    def __init__ ( self, M0, easy=None, difficult=None ):
        self.hf0 = M0.hf0
        self.w   = M0.w.copy()
        self.nu  = M0.nu
        self.easy = easy
        self.difficult = difficult
    def __call__ ( self, fitted_model ):
        # We want to use the fitted model
        self.hf0 = fitted_model.hf0
        self.w = fitted_model.w
        self.nu = fitted_model.nu
        self.pi = fitted_model.pi

        X_ = fitted_model.X.copy ()
        for j in fitted_model.applythreshold:
            X_[:,j] = threshold.u_v ( X_[:,j], self.nu )
        self.history_evaluation ( X_, fitted_model.r, (0.75,.55) )

        out = [Likelihood_and_Slope_Collector ( fitted_model )] #0,1,2
        out.append ( [self.vdifficult] ) # 3
        out.append ( [self.veasy] )      # 4
        out.append ( [self.phist] )      # 5
        out.append ( [self.pstim] )      # 6
        out.append ( [self.pSH]   )      # 7
        out.append ( [self.peasy] )      # 8
        out.append ( fitted_model.pi ) # 9,10,11
        # out.append ( fitted_model.w[1:self.hf0] ) # 12,...,12+hf0
        out.append ( self.get_thres ( fitted_model.w[:self.hf0], fitted_model.nu, fitted_model.pi, .75 ) ) # 12,...,12+hf0
        out.append ( self.get_thres ( fitted_model.w[:self.hf0], fitted_model.nu, fitted_model.pi, .85 ) ) # 12+hf0,...,12+2*hf0
        out.append ( fitted_model.nu ) # 
        out.append ( fitted_model.w ) # 
        out.append ( np.ravel ( self.stimuli ) )
        out.append ( np.ravel ( self.conditions ) )
        out.append ( np.ravel ( self.variance_explained ) )
        return np.concatenate ( out )

    def get_thres ( self, w, nu, pi, p=.75 ):
        """Determine the p-Threshold on performance correct"""
        out = []
        for i in xrange ( 1, w.shape[0] ):
            def e ( x ):
                p1 = pi[1]+pi[2]*model.logistic ( w[0]+w[i]*threshold.u_v ( x, nu ) )
                p2 = 1-(pi[1]+pi[2]*model.logistic ( w[0]+w[i]*threshold.u_v ( -x, nu ) ))
                return abs(.5*p1+.5*p2-p)
            out.append ( float(opt.fmin ( e, 10., disp=0 )) )
        return np.array(out)

    def history_evaluation ( self, X, r, p ):
        """Determine variance explained by history

        :Parameters:
            *X*
                design matrix
            *r*
                responses
            *p*
                probability correct that is considered the border between
                easy and difficult
        """
        if self.easy is None or self.difficult is None:
            difficult,easy = model.performance_filter ( r, X, p[0], p[1], hf0=self.hf0 )
        else:
            difficult,easy = self.difficult,self.easy
        # easy = np.logical_not ( difficult )

        current_stimulus = np.dot ( X[:,1:self.hf0], self.w[1:self.hf0] )
        history_features = np.dot ( X[:,self.hf0:], self.w[self.hf0:] )
        decision_signal  = current_stimulus + history_features
        difficult_var = np.var ( current_stimulus[difficult] ) + \
                np.var ( history_features[difficult] )
        easy_var = np.var ( current_stimulus[easy] ) + \
                np.var ( history_features[easy] )

        self.vdifficult = np.var ( history_features[difficult] ) / \
                difficult_var
        self.veasy      = np.var ( history_features[easy] ) / easy_var
        self.vstimulus  = np.var ( current_stimulus[easy] ) / easy_var

        S = []
        V = []
        C = []

        for condition in xrange ( 1, self.hf0 ):
            stimuli = np.unique ( abs(X[:,condition]) )
            S_ = []
            V_ = []
            for s in stimuli:
                if abs(s)<1e-10:
                    continue
                i = abs(X[:,condition]) == s
                S_.append ( s )
                V_.append ( np.var(history_features[i])/np.var(history_features[i]+current_stimulus[i]) )
            S.append ( S_ )
            V.append ( V_ )
            C.append ( [condition]*len(S_) )
        self.stimuli            = np.concatenate ( S )
        self.conditions         = np.concatenate ( C )
        self.variance_explained = np.concatenate ( V )

        pred_hist = self.pi[1]+self.pi[2]*model.logistic ( history_features[difficult] + self.w[0] )
        pred_stim = self.pi[1]+self.pi[2]*model.logistic ( current_stimulus[difficult] + self.w[0] )
        pred_SH   = self.pi[1]+self.pi[2]*model.logistic ( decision_signal[difficult]  + self.w[0] )
        pred_easy = self.pi[1]+self.pi[2]*model.logistic ( history_features[easy] + self.w[0] )
        self.phist = np.mean ( (pred_hist>0.5) == r[difficult] )
        self.pstim = np.mean ( (pred_stim>0.5) == r[difficult] )
        self.pSH   = np.mean ( (pred_SH  >0.5) == r[difficult] )
        self.peasy = np.mean ( (pred_easy>0.5) == r[easy] )

#################### The real monte carlo stuff ####################

def mcsimulation ( generator, collector, nsamples, w0, p0, nu0=0., verbose=False, hf0=2, applythreshold=[] ): 
    """run a monte carlo simulation

    :Parameters:
        *generator*
            a callable such that generator() will generate a new pair of
            responses and design matrix
        *collector*
            a callable such that collector(fitted_model) will be everything
            that is stored from the fitted dataset
        *nsamples*
            how many samples are to be generated
        *w0*
            starting values for the optimizations
        *verbose*
            if true, print status bar

    :Returns:
        a list of the collector outputs

    :Example:
    >>> np.random.seed(0)
    >>> generator = lambda: ((np.random.rand(10)>0.5).astype('d'),np.c_[np.ones(10),np.random.uniform(-1,1,size=10)])
    >>> collector = lambda M: (M.w,M.pi)
    >>> sim_results = mcsimulation ( generator, collector, 3, [10,10] )
    >>> len(sim_results)
    3
    >>> sim_results[0]
    (array([ 0.89742561,  1.36421569]), array([ 0.04270847,  0.03216249,  0.92512904]))
    """
    if verbose:
        status = np.arange ( 0, nsamples, nsamples / 10 )
    else:
        status = [-1]
    collected_output = []
    for i in xrange ( nsamples ):
        r, X = generator ()
        # print(np.shape(r))
        # print(np.shape(X))
        
        M = model.history_model ( r, X, w0=w0, p0=p0, nu0=nu0, hf0=hf0,
                applythreshold=applythreshold, emiter=80 )
        collected_output.append ( collector(M) )
        w0 = M.w
        p0 = M.pi
        nu0 = M.nu
        if i in status:
            if verbose:
                sys.stderr.write ( "\r%d%%" % (100*(float(i)/nsamples), ) )
                sys.stderr.flush()
    if verbose:
        sys.stderr.write ( "\rDone\n" )
        sys.stderr.flush()
    return collected_output

#################### Take some special collected data and derive something meaningful from them ####################

def historytest ( simulated_likelihoods, nbins=20 ):
    """perform a likelihood based test for history dependency

    :Parameters:
        *simulated_likelihoods*
            an array of loglikelihoods on permuted data
        *nbins*
            number of histogram bins

    :Returns:
        (h,b)
            a histogram of the loglikelihoods
        C95
            critical value from the 95th percentile of the likelihood distribution
    """
    hist = np.histogram ( simulated_likelihoods, bins=nbins )
    C95 = prctile ( simulated_likelihoods, 0.95 )
    return hist,C95

def history_kernel_ci ( stimulus_kernels, response_kernels, point_stimulus, point_response ):
    """determine confidence intervals for history kernels using bootstrap

    :Parameters:
        *stimulus_kernels*
            an array with shape (nsamples,nlags) in which each row corresponds
            to one simulated stimulus kernel
        *response_kernels*
            an array with shape (nsamples,nlags) in which each row corresponds
            to one simulated response kernel

    :Returns:
        ci_stimulus,ci_response,ci_correct,ci_error
            each of these corresponds to a (2,nlags) array of upper and lower
            confidence limits

    :Example:
    >>> stimulus_kernels = np.array([np.mgrid[0:1:1001j] for i in xrange (3)]).T
    >>> response_kernels = np.array([np.mgrid[0:1:1001j] for i in xrange (3)]).T
    >>> ci_k = history_kernel_ci ( stimulus_kernels, response_kernels )
    >>> ci_k[0]
    array([[ 0.025,  0.025,  0.025],
           [ 0.975,  0.975,  0.975]])
    >>> ci_k[2]
    array([[ 0.05,  0.05,  0.05],
           [ 1.95,  1.95,  1.95]])
    """
    stimulus_kernels = np.array(stimulus_kernels)
    response_kernels = np.array(response_kernels)
    nsamples,nlags = stimulus_kernels.shape
    assert response_kernels.shape == (nsamples,nlags)

    ci_stimulus,ci_response = np.zeros ( (2,nlags), 'd' ), np.zeros ( (2,nlags), 'd' )
    ci_correct,ci_error     = np.zeros ( (2,nlags), 'd' ), np.zeros ( (2,nlags), 'd' )

    correct_kernels =  stimulus_kernels+response_kernels
    error_kernels   = -stimulus_kernels+response_kernels
    point_correct   =  point_stimulus+point_response
    point_error     = -point_stimulus+point_response

    ci_method = 'plugin'

    for lag in xrange ( nlags ):
        ci_stimulus[:,lag] = bootstrap_ci ( point_stimulus[lag],stimulus_kernels[:,lag], method=ci_method )
        ci_response[:,lag] = bootstrap_ci ( point_response[lag],response_kernels[:,lag], method=ci_method )
        ci_correct [:,lag] = bootstrap_ci ( point_correct[lag], correct_kernels[:,lag], method=ci_method )
        ci_error   [:,lag] = bootstrap_ci ( point_error[lag],   error_kernels[:,lag], method=ci_method )

    return ci_stimulus,ci_response,ci_correct,ci_error

def bootstrap_ci ( xhat, x_, method='plugin' ):
    if method == 'plugin':
        ci = prctile ( x_, (.025, .975) )
    elif method == 'se':
        se = np.std ( x_ )
        ci = xhat - 2*se, xhat + 2*se
    elif method == 'pivot':
        prc = prctile ( x_, (.975,.025) )
        ci = 2*xhat-prc[0],2*xhat-prc[1]
    else:
        raise NotImplementedError
    return ci

#################### Helper functions ####################

def prctile ( x, fractions ):
    """percentiles of an array

    :Parameters:
        *x*
            an array containing the samples
        *fractions*
            either a float or a sequence of floats marking the desired fractions

    :Examples:
    >>> x = [1,2,3,4,5]
    >>> prctile ( x, .9 )
    5
    >>> prctile ( x, .5 )
    3
    >>> prctile ( x, (.5,.9) )
    array([3, 5])
    """
    x_sorted = np.sort ( x )
    # if not getattr ( fractions, '__iter__', False ):
    #     fractions = [fractions]
    indices = np.ceil(np.array(fractions)*(len(x)-1)).astype('i')
    out = x_sorted[indices]
    return out

def aic ( independent_model, history_model ):
    """Determine critical likelihood value based on Akaike's information criterion

    :Parameters:
        *independent_model*
            model with no history terms
        *history_model*
            model with history terms
    """
    prmdiff = len(history_model.w)-len(independent_model.w)
    l0 = independent_model.loglikelihood
    return l0 + prmdiff

if __name__ == "__main__":
    import doctest
    doctest.testmod()
