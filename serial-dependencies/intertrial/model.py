#!/usr/bin/env python

import numpy as np
import threshold,glm
import sys

__doc__ = """This function defines the actual history model and allows to fit this model to data.

Copyright (C) 2014 Ingo Fruend

This code reproduces the analyses in the paper

    Fruend, Wichmann, Macke (2014): Quantifying the effect of inter-trial dependence on perceptual decisions. J Vis, 14(7): 9.


    Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

    If you use the Software for your own research, cite the paper.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""
__all__ = ["history_model"]

class history_model ( object ):
    def __init__ ( self, r, X, **kwargs ):
        """History dependent model

        :Parameters:
            *r*
                responses
            *X*
                design matrix
            *hf0*
                index of the first history features
            *nafc*
                number of stimulus alternatives presented.
            *applythreshold*
                list of features to which the thresholding nonlinearity should be applied.
            *verbose*
                show starting values if this is true

        :Starting values and priors:
            *pprior*
                'prior' for the mixing coefficients. These values can be interpreted in terms of
                the parameter vector for a Dirichlet-prior on the mixing coefficients. If ``pprior``
                is a vector b, then the mixing coeffients have a prior Dir(b+1).
            *w0*
                starting values for the glm weights w. If this is shorter than the actual vector to
                be fitted, zeros will be appended at the end.
            *p0*
                starting values for the mixing coefficients pi. The first value corresponds to
                the probability to guess a 1 response, the second value corresponds to the probability
                to guess a 0 response and the third value corresponds to the probability to respond
                based on the stimulus.
            *nu0*
                starting value for the threshold. The threshold will only be optimized if ``applythreshold``
                is not empty

        :Parameters for the algorithm:
            *emiter*
                number of EM iterations to be performed
            *miniter*
                minimum number of EM iterations to be performed (to make sure that the model
                has time to move away from initial shallow regions that are due to bad starting
                values)
            *nuiter*
                number of iterations in optimizing the threshold. This refers to a single EM step.
            *nustop*
                stopping criterion for nu on each EM step
            *glmiter*
                number of IRLS iterations when fitting the latent generalized linear model.
            *glmstop*
                stopping criterion for the glm iterations

        :Example:
        Not a particularly good one in fact
        >>> r = np.array([1, 1, 1, 1, 0, 1, 0, 1, 1, 0], dtype='d')
        >>> X = np.c_[np.ones(len(r)),np.array([0.58345008,0.05778984,0.13608912,0.85119328,-0.85792788,-0.8257414,-0.95956321,0.66523969,0.5563135,0.7400243])]
        >>> M = history_model ( r, X )
        """
        self.hf0  = kwargs.setdefault ( "hf0",  2 )
        self.nafc = kwargs.setdefault ( "nafc", 2 )

        self.nuiter  = kwargs.setdefault ( 'nuiter',  10 )
        self.glmiter = kwargs.setdefault ( 'glmiter', 20 )
        self.emiter  = kwargs.setdefault ( 'emiter',  20 )
        self.miniter = kwargs.setdefault ( 'miniter', 10 )

        self.nustop  = kwargs.setdefault ( 'nustop',  1e-6 )
        self.glmstop = kwargs.setdefault ( 'glmstop', 1e-6 )
        self.emstop  = kwargs.setdefault ( 'emstop',  1e-6 )
        self.emabs   = kwargs.setdefault ( 'emabs',   1e-4 )
        self.storeopt= kwargs.setdefault ( 'storeopt',False )
        self.pprior  = np.array(kwargs.setdefault ( 'pprior',  [0]*3 ),'d')
        self.lmprior = kwargs.setdefault ( 'lm', 0.1 )

        self.applythreshold = kwargs.setdefault ( 'applythreshold', [] )

        self.verbose = kwargs.setdefault ( 'verbose', False )
        self.verbose = False # no output

        # # This would allow setting the starting parameters
        w0, p0, nu0 = history_model.__heuristics_for_start ( r, X, self.hf0, self.nafc )
        self.w = np.zeros ( X.shape[1], 'd' )
        w  = kwargs.setdefault ( 'w0', w0 )
        
        # in the modulation and 7 independent lags model, this goes wrong
        try:
            self.w[:len(w0)] = w0
        except:
            from IPython import embed
            embed()    
            
        self.w0 = self.w.copy()
        self.pi = kwargs.setdefault ( 'p0', p0 )
        self.nu = kwargs.setdefault ( 'nu0', nu0 )

        if len ( self.applythreshold ) == 0:
            self.nu = 0.
        if self.storeopt:
            self.opt = []

        self.w, self.pi, self.q, self.loglikelihood, self.nu = self.__em ( X, r )

        # self.history_evaluation ( X, r, p=0.75 )
        self.X = X
        self.r = r

    def __em ( self, X, r ):
        """Optimize parameters using expectation maximation

        :Parameters:
            *X*
                design matrix
            *r*
                response vector
        """
        w  = self.w
        p  = self.pi
        nu = self.nu
        # w = np.array ( [.01,.8,-.8,.5] )

        if self.verbose:
            print "Starting values:"
            print "nu:",nu
            print "w: ",w
            print "p: ",p

        X_ = X.copy ()
        for j in self.applythreshold:
            X_[:,j] = threshold.u_v ( X[:,j], nu )

        # Expectation
        gwx = history_model.__combine_features ( X_, w )
        q = history_model.__determine_single_trial_lapses ( r, gwx, p )
        p = history_model.__optimize_p ( q, self.pprior )

        l = history_model.__likelihood ( gwx, r, p )

        for i in xrange ( self.emiter ):
            # Maximization
            nu_ = nu
            nu = threshold.optimize_nu ( X, r, q[:,-1], w, nu, self.applythreshold, niter=self.nuiter, stop=self.nustop )
            if np.isnan ( nu ):
                sys.exit ( 2 )
            for j in self.applythreshold:
                X_[:,j] = threshold.u_v ( X[:,j], nu )
            # if len ( self.applythreshold ) == 1:
            #     nu = threshold.optimize_nu ( X, r, q[:,-1], w, nu, niter=self.nuiter, stop=self.nustop )
            #     if np.isnan(nu):
            #         sys.exit ( 2 )
            #     for j in self.applythreshold:
            #         X_[:,j] = threshold.u_v ( X[:,j], nu )
            # elif len ( self.applythreshold ) > 1:
            #     raise NotImplementedError, "optimization of the threshold nu is not implemented for this case"

            w_ = w
            p_ = p
            w = glm.optimize_w ( X_, r, q[:,-1], w, niter=self.glmiter, stop=self.glmstop, lm=self.lmprior )
            p = history_model.__optimize_p ( q, self.pprior )
            if np.isnan ( p ).any():
                print i,p
                sys.exit ( 1 )

            # Expectation
            gwx = history_model.__combine_features ( X_, w )
            q   = history_model.__determine_single_trial_lapses ( r, gwx, p )

            l_ = history_model.__likelihood ( gwx, r, p )

            # Stop?
            rel_e = np.abs ( (l_-l)/l )
            abs_e = max ( max ( np.abs(w-w_).max(),np.abs(p-p_).max() ), abs(nu-nu_) )
            if i>self.miniter:
                if rel_e < self.emstop and abs_e < self.emabs:
                    if self.verbose:
                        sys.stderr.write ( "Converged after %d iterations\n  relative error: %g\n  absolute error: %g\n" % (i,rel_e,abs_e) )
                    break
            if self.storeopt:
                self.opt.append ( [l_,rel_e,abs_e]+[ pp for pp in p]+[ ww for ww in w ]+[float(nu)] )

            if self.verbose:
                print l_,rel_e,abs_e
            l = l_
        else:
            if self.verbose:
                sys.stderr.write ( "No convergence after %d iterations\nrelative error: %g\n" % (i,rel_e ) )

        return w, p, q, l, nu

    @staticmethod
    def __heuristics_for_start ( r, X, hf0, nafc ):
        """Determine starting values for optimization based on some heuristics

        :Parameters:
            *r*
                vector of responses
            *X*
                design matrix
            *hf0*
                index of the column with the first history feature
            *nafc*
                number of forced choice alternatives
        """
        w0 = np.zeros ( X.shape[1], 'd' )
        for j in xrange ( 1, hf0 ):
            smax = X[:,j].max()
            if nafc==2:
                gt0 = X[:,j]>0
                nu0 = (X[gt0,j]*r[gt0]).mean()
                m = X[:,0].mean()
                smax -= nu0
                smin = -smax
            else:
                nu0 = 0.
                m = np.unique ( X[:,j] ).mean ()
                smin = X[:,j].min()
            w = 4*np.log(5)/(smax-smin)
            w0[j] = w
        w0[0] = -m*w
        p0 = [.01,.01,.98]
        # p0 = [1e-5,1e-5,1.-2e-5]
        return w0,p0,np.sqrt(nu0)

    @staticmethod
    def __combine_features ( X_thres, w ):
        """Combine features linearly and apply logistic

        :Parameters:
            *X_thres*
                design matrix (after potential application of the threshold)
            *w*
                feature weights
        """
        eta = np.dot ( X_thres, w )
        return logistic ( eta )

    @staticmethod
    def __likelihood ( gwx, r, p ):
        """Determine the log likelihood of the data r

        :Parameters:
            *gwx*
                response probabilities if the observer would look at the stimulus
                This is a linear combination of the features mapped through a logistic
                function.
            *r*
                response vector
            *q*
                prior state probabilities
        """
        # prob0 = p[0] + p[2]*(1-gwx)
        prob1 = p[1] + p[2]*gwx
        return np.sum ( r*np.log(prob1) + (1-r)*np.log1p(-prob1) )

    @staticmethod
    def __optimize_p ( q, prior ):
        """Optimize global state probabilities

        :Parameters:
            *q*
                previous single trial state probabilities
            *prior*
                parameters of a possible dirichlet prior on the state probabilities (1 is added to each of these)
        """
        # Closed form solution is available
        pnew = q.sum(0)
        pnew += prior
        return pnew / q.shape[0]

    @staticmethod
    def __determine_single_trial_lapses ( r, gwx, p ):
        """Determine the single trial posterior distribution q of lapses

        :Parameters:
            *r*
                single trial binary responses
            *gwx*
                response probabilities if the observer would look at the stimulus
                This is a linear combination of the features mapped through a logistic
                function.
            *p*
                current global state probabilities

        :Example:
        >>> r = [0,1,1]
        >>> gwx = [.7,.5,.7]
        >>> p = [.1,.1,.8]
        >>> q = history_model._history_model__determine_single_trial_lapses ( r, gwx, p )
        >>> q[0,1] == 0 and q[1,0] == 0.0 and q[2,0] == 0.0
        True
        >>> print round ( q[0,0], 2 ), round ( q[0,2], 2 )
        0.29 0.71
        >>> print round ( q[1,1], 2 ), round ( q[1,2], 2 )
        0.2 0.8
        >>> print round ( q[2,1], 2 ), round ( q[2,2], 2 )
        0.15 0.85
        """
        r = np.array(r)
        p = np.array(p)
        gwx = np.array(gwx)
        # P(q|r) = P(r|q) * P(q) / Z = P(r|q) * p[q] / Z
        q0 = (1.-r)*p[0]
        q1 = r*p[1]
        q2 = r*gwx + (1-r)*(1-gwx)
        q2 *= p[2]
        q = np.c_[q0,q1,q2]
        q /= np.reshape ( q.sum(1),(-1,1) )
        return q

    @staticmethod
    def evaluate ( X, w, p ):
        """Evaluate the model on the design matrix X

        :Parameters:
            *X*
                Design matrix with one trial per row and one feature per column
            *w*
                feature weights
            *p*
                state probabilities
        :Example:
        >>> X = np.ones ( (4,2), 'd' )
        >>> X[:,1] = np.arange ( 4 )
        >>> w = [.5,.5]
        >>> p = [.05,.05,.9]
        >>> history_model.evaluate ( X, w, p )
        array([ 0.6102134 ,  0.70795272,  0.78581703,  0.84271737])
        """
        assert np.abs(np.sum(p)-1)<1e-8
        Xw = np.dot ( X,w )
        return p[1] + p[2]*logistic ( Xw )

    def __getitem__ ( self, key ):
        return self.__dict__[key]

    def pcorrect ( self, x, pleft=0.5, pright=0.5, ind=[0,1] ):
        """Get probability of a correct response rather than probability of left/right response
        
        :Parameters:
            *x*         stimulus intensities
            *pleft*     probability that the stimulus is on the left (or number of stimuli on the left)
            *pright*    probability that the stimulus is on the right
        """
        Z = pleft+pright
        if Z>1:
            pleft /= Z
            pright /= Z
        s = threshold.u_v(x,self.nu)
        psi_p = self.pi[1]+self.pi[2]*logistic (    self.w[ind[0]]+self.w[ind[1]]*s )
        psi_m = self.pi[0]+self.pi[2]*(1-logistic ( self.w[ind[0]]-self.w[ind[1]]*s ))
        return pleft*psi_m + pright*psi_p,s

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
        difficult = performance_filter ( r, X, p, hf0=self.hf0 )
        easy = np.logical_not ( difficult )

        X_ = X.copy ()
        for j in self.applythreshold:
            X_[:,j] = threshold.u_v ( X[:,j], self.nu )

        current_stimulus = np.dot ( X_[:,1:self.hf0], self.w[1:self.hf0] )
        history_features = np.dot ( X_[:,self.hf0:], self.w[self.hf0:] )

        self.vdifficult = np.var ( history_features[difficult] ) / \
                (np.var (history_features[difficult]) + np.var(current_stimulus[difficult]))

        self.veasy      = np.var ( history_features[easy] ) / \
                (np.var (history_features[easy]) + np.var(current_stimulus[easy]))

        S = []
        V = []
        C = []

        for condition in xrange ( 1, self.hf0 ):
            stimuli = np.unique ( abs(X_[:,condition]) )
            S_ = []
            V_ = []
            for s in stimuli:
                if abs(s)<1e-10:
                    continue
                i = abs(X_[:,condition]) == s
                S_.append ( s )
                V_.append ( np.var(history_features[i])/(np.var(history_features[i])+np.var(current_stimulus[i])) )
            S.append ( S_ )
            V.append ( V_ )
            C.append ( [condition]*len(S_) )
        self.stimuli            = np.concatenate ( S )
        self.conditions         = np.concatenate ( C )
        self.variance_explained = np.concatenate ( V )

        # c = -self.w[0]
        c = 0.
        self.phist = np.mean (
                (history_features[difficult]>c) == r[difficult] )
        self.pstim = np.mean (
                (current_stimulus[difficult]>c) == r[difficult] )
        self.pSH   = np.mean ( (
            (history_features[difficult]+current_stimulus[difficult])>c )
            == r )
        self.peasy = np.mean ( (history_features[easy]>c) == r[easy] )

        return np.var ( current_stimulus ), np.var ( history_features ),(S,V)

def logistic ( x ):
    """Logistic function"""
    return 1./( 1.+np.exp(-x) )

def performance_filter ( r, X, p1=0.75, p0=0.55, hf0=2 ):
    """Filter data based on performance

    :Parameters:
        *r*
            response vector
        *X*
            design matrix
        *p1*
            upper performance level
        *p0*
            lower performance level
        *hf0*
            index of first history feature

    :Returns:
        difficult,easy
            indices into X and r such that r[difficult],X[difficult] are
            stimulus levels for which performance between p0 and p1 is
            expected and r[easy],X[easy] are stimulus levels for which
            performance better than p1 is expected.
    """
    stimulus = X[:,1:hf0].sum(1)
    target   = stimulus >= 0
    correct  = target==r
    stim_levels = np.unique ( abs ( stimulus ) )

    difficult = np.zeros ( X.shape[0], 'bool' )
    easy      = np.zeros ( X.shape[0], 'bool' )

    for c in xrange ( 1, hf0 ):
        for s in stim_levels:
            index = abs(X[:,c])==s
            pcorrect = np.mean ( correct[index] )
            if pcorrect > p1:
                easy = np.logical_or ( easy, index )
            elif pcorrect > p0:
                difficult = np.logical_or ( difficult, index )

    return difficult,easy

if __name__ == "__main__":
    import doctest
    doctest.testmod ()
