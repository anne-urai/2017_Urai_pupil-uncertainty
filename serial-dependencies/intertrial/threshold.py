#!/usr/bin/env python

import numpy as np
import scipy.optimize as opt
import sys,cPickle
import scipy.weave

__doc__ = """Definitions for the threshold nonlinearity


Copyright (C) 2014 Ingo Fruend

This code reproduces the analyses in the paper

    Fruend, Wichmann, Macke (2014): Quantifying the effect of inter-trial dependence on perceptual decisions. J Vis, 14(7): 9.


    Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

    If you use the Software for your own research, cite the paper.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

# These two are fixed in the paper.
thres = True # Use threshold nonlinearity (rather than power)
kappa = 4.   # Use threshold exponent of 4 (relatively sharp threshold)

weave_code_uxw = r"""
int i,j;
bool t;
double eta,u,exnu,xij,wj;

for (i=0; i<X_array->dimensions[0]; i++) {
    eta = 0;
    for ( j=0; j<X_array->dimensions[1]; j++) {

        xij = *((double*)(X_array->data+i*X_array->strides[0]+j*X_array->strides[1]));
        wj  = *((double*)(w_array->data+j*w_array->strides[0]));
        t   = *((int*)(thres_array->data+j*thres_array->strides[0]));

        if ( t ) {
            exnu = exp(xij-nu2);
            u = log(1+exnu*exnu*exnu*exnu);
            exnu = exp(-xij-nu2);
            u -= log(1+exnu*exnu*exnu*exnu);
            eta += wj*u/4.;
        } else {
            eta += wj*xij;
        }
    }

    *(double*)(out_array->data+i*out_array->strides[0]) = eta;
}
"""

def g ( x ):
    """The logistic function"""
    return 1./(1+np.exp(-x))

def dg ( gxw ):
    """derivative of the psychometric function"""
    return gxw*(1-gxw)

def ddg ( gxw ):
    """2nd derivative of the psychometric function"""
    return gxw*(1-gxw)*(1-2*gxw)

if thres:
    def u_v ( x, nu ):
        """The soft threshold function"""
        return (np.log(1+np.exp(x-nu**2)**kappa)-np.log(1+np.exp(-x-nu**2)**kappa))/kappa
    def dudnu ( x, nu ):
        """Derivative of u_v with respect to nu"""
        return 2*nu*(g(-kappa*(x-nu**2)) - g( kappa*(x+nu**2)))
    def ddudnu ( x, nu ):
        """2nd derivative of u_v with respect to nu"""
        sg = g( kappa*(x+nu**2))
        sg_= g(-kappa*(x-nu**2))
        return 4 * kappa * nu**2 * (sg_*(1-sg_) - sg*(1-sg)) + 2 * (sg_-sg)
else:
    def u_v ( x, nu ):
        """An accelerating transducer function"""
        return np.sign ( x ) * np.abs ( x )**nu

    def dudnu ( x, nu ):
        """Derivative of u_v with respect to nu"""
        return np.sign ( x ) * np.abs ( x )**nu * np.log ( np.abs(x) )

    def ddudnu ( x, nu ):
        """2nd derivative of u_v with respect to nu"""
        return np.sign ( x ) * np.abs ( x )**nu * np.log ( np.abs(x) )**2

# Evaluation of the whole model can be accellerated using compiled code
def psi_weave ( X, w, nu, applythreshold ):
    """Evaluation of the whole model using compiled code"""
    assert w.shape[0] == X.shape[1]
    thres = np.zeros ( w.shape[0], np.int32 )
    thres[applythreshold] = 1
    out = np.zeros ( X.shape[0], np.float64 )
    nu2 = float(nu*nu)
    scipy.weave.inline ( weave_code_uxw, ['X','w','nu2','thres','out'] )
    return g(out)

def psi_py ( X, w, nu, applythreshold ):
    """Evaluation of the whole model using only python"""
    eta = w[0]*X[:,0]
    for j in xrange ( 1, X.shape[1] ):
        if j in applythreshold:
            eta += w[j]*u_v ( X[:,j], nu )
        else:
            eta += w[j]*X[:,j]
    return g ( eta )

# psi is clearly a bottleneck in profiling the code.
# it is highly recommended to use the weave based evaluation of psi
psi = psi_py

def dpsidnu ( gxw, w, dudnu_ ):
    """derivative of psi with respect to nu"""
    return dg(gxw)*w[1]*dudnu_

def ddpsidnu ( gxw, w, dudnu_, ddudnu_ ):
    """2nd derivative of psi with respect to nu"""
    return ddg(gxw)*w[1]**2*dudnu_**2 + dg(gxw)*w[1]*ddudnu_

def L ( gxw, r, q ):
    """log likelihood"""
    return np.sum ( q*(r*np.log(gxw) + (1-r)*np.log(1-gxw)) )

def dL ( gxw, r, q, dpsidnu_ ):
    """derivative of log likelihood"""
    return np.sum( q*dpsidnu_ * ( np.clip(np.where ( r==0, 0, r/gxw ) - np.where ( r==1, 0, (1-r)/(1-gxw) ), 0, 1e5 ) ) )

def ddL ( gxw, r, q, dpsidnu_, ddpsidnu_ ):
    """second derivative of log likelihood"""
    return np.sum(q*(
            ddpsidnu_ * (np.clip(np.where ( r==0, 0, r/gxw ) - np.where ( r==1, 0, (1-r)/(1-gxw)), 0, 1e5)) \
                    - dpsidnu_**2 * ( np.clip(np.where ( r==0, 0, r/gxw**2 ) + np.where ( r==1, 0, (1-r)/(1-gxw)**2 ), 0, 1e10 ) ) ) )

def optimize_nu ( X, y, q, w, nu0, applythreshold=[1], niter=10, stop=1e-5, longoutput=False ):
    """optimize nu"""
    def error ( nu ):
        gxw = psi ( X, w, nu, applythreshold )
        return - L ( gxw, y, q )
    nu_out = opt.fmin ( error, float ( nu0 ), disp=0 )
    return nu_out

