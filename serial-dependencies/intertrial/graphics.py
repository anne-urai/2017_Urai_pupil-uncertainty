#!/usr/bin/env python

import pylab as pl
import numpy as np
from matplotlib.patches import Arc
from matplotlib.ticker import MaxNLocator
from matplotlib import gridspec
from scipy import stats

__doc__ = """Graphics commands that allow generating all plots in the canonical figure


Copyright (C) 2014 Ingo Fruend

This code reproduces the analyses in the paper

    Fruend, Wichmann, Macke (2014): Quantifying the effect of inter-trial dependence on perceptual decisions. J Vis, 14(7): 9.


    Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

    If you use the Software for your own research, cite the paper.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

# A number of settings for the homogeneous appearance of all plots
histogram_color = np.array ( [189,189,189] ) / 255. # gray 74
observed_color  = np.array ( [  0,  0, 88] ) / 255. # gray 12
C95_color       = np.array ( [150,150,150] ) / 255.
Caic_color      = np.array ( [139, 30, 30] ) / 255.
stimulus_color  = np.array ( [ 25, 25,112] ) / 255.
response_color  = np.array ( [255,215,  0] ) / 255.
correct_color   = np.array ( [ 34,139, 34] ) / 255.
incorrect_color = np.array ( [240,255,255] ) / 255.
linewidths      = 1

class canonical_axes:
    def __init__ ( self, fig, bbox=(0,0,1,1) ):
        """Determine the general layout of the axes"""
        x0,y0,w,h = bbox
        gs = gridspec.GridSpec(3, 3)
        self.pmf               = fig.add_subplot(gs[0,0])
        self.uv                = fig.add_subplot(gs[0,1])
        self.likeli            = fig.add_subplot(gs[0,2])
        self.history_rz        = fig.add_subplot(gs[1,0])
        self.history_perf      = fig.add_subplot(gs[1,1])
        self.slopes            = fig.add_subplot(gs[1,2])
        self.history_rzP       = fig.add_subplot(gs[2,0])
        self.history_perfP     = fig.add_subplot(gs[2,1])
        self.history_pupil     = fig.add_subplot(gs[2,2])
        

def prepare_axes ( ax, haveon=("bottom","left" ) ):
    """Prepare an axes object to look nicer than standard matplotlib

    :Parameters:
        *ax* :
            axes object that should be prepared
        *haveon* :
            axes that should be shown

    :Return:
        the prepared axes object
    """
    if getattr ( ax, 'spines', False ):
        for loc,spine in ax.spines.iteritems():
            if loc in haveon:
                spine.set_position ( ("outward",3) )
            else:
                spine.set_color ( "none" )
    else:
        warnings.warn ( spineswarning, DeprecationWarning )

    if "bottom" in haveon:
        ax.xaxis.set_ticks_position ( "bottom" )
    elif "top" in haveon:
        ax.xaxis.set_ticks_position ( "top" )
    else:
        ax.xaxis.set_ticks_position ( "none" )
        ax.xaxis.set_ticklabels ( "" )
    if "left" in haveon:
        ax.yaxis.set_ticks_position ( "left" )
    elif "right" in haveon:
        ax.yaxis.set_ticks_position ( "right" )
    else:
        ax.yaxis.set_ticks_position ( "none" )
        ax.yaxis.set_ticklabels ( "" )

    return ax

def label_axes ( title=None, xlabel=None, ylabel=None, legend=None, ax=None, nxticks=None, nyticks=None ):
    """Set all labels for an axes

    :Parameters:
        *title*
            title string
        *xlabel*
            x-axis label
        *ylabel*
            y-axis label
        *legend*
            position of the legend (the legend is created from plot labels!)
        *ax*
            pylab axes (if None, this will by pylab.gca())

    :Example:
    >>> l = pl.plot ( [1,3,2], label='test' )
    >>> label_axes ( title='test', xlabel=r'$x$', ylabel=r'$f(x)$', legend='lower right' )
    >>> pl.savefig ( 'test/label_axes.png' ); pl.close()
    """
    if ax is None:
        ax = pl.gca ()

    if not title is None:
        ax.set_title ( title )
    if not xlabel is None:
        ax.set_xlabel ( xlabel )
    if not ylabel is None:
        ax.set_ylabel ( ylabel )
    if not legend is None:
        ax.legend ( loc=legend, numpoints=1,
            handletextpad=.3, columnspacing=.15, borderpad=.2,
            handlelength=1, labelspacing=.1
            )

def montecarlo_test ( observed_value, hist, C95, Caic=None, ax=None, labeling="likelihood" ):
    """Graphical comparison of an observed value to a histogram and a critical value

    :Parameters:
        *observed_value*
            the observed value
        *hist*
            histogram of simulated values
        *C95*
            Criterion from the 95th percentile of the histogram distribution
        *Caic*
            Criterion based on Akaike's information criterion
        *ax*
            pylab.axes object in which the plot should go
        *labeling*
            either 'likelihood' or 'slope'

    :Example:
    >>> x = pl.randn(100)
    >>> obs = 3.
    >>> hist = np.histogram ( x )
    >>> C95 = pl.prctile ( x, 95 )
    >>> Caic = 2.5
    >>> montecarlo_test ( obs, hist, C95, Caic )
    >>> pl.savefig ( 'test/montecarlo_test.png' ); pl.close()
    """
    if ax is None:
        ax = pl.gca()
    ax = prepare_axes ( ax, haveon=("bottom",) )

    if not hist is None:
        h,b = hist
        xti = [np.floor(np.min(b)),np.ceil(np.max(b))]
        ax.bar ( b[:-1], h, pl.diff(b), edgecolor=histogram_color, facecolor=histogram_color )
        yrange = ax.get_ylim ()
    else:
        yrange = (np.min([observed_value,Caic])-1,np.max([observed_value,Caic])+1)

    if not labeling in ["likelihood","slope","other"]:
        raise ValueError, "labeling should be either 'likelihood' or 'slope'"

    if observed_value>C95:
        ax.plot ( [observed_value]*2, (yrange[0],yrange[0]+0.85*(yrange[1]-yrange[0])),
            color=observed_color, linewidth=2,
            label=r"$\ell_\mathrm{obs}$" if labeling=="likelihood" else r"$\psi'(s)$" )
        ax.plot ( [observed_value], [yrange[0]+0.95*(yrange[1]-yrange[0])], '*', color=observed_color )
        # ax.text ( observed_value, yrange[0]+0.95*(yrange[1]-yrange[0]), "*",
        #     horizontalalignment='center', verticalalignment='center', fontsize=12 )
    else:
        ax.plot ( [observed_value]*2, yrange,
            color=observed_color, linewidth=2,
            label=r"$\ell_\mathrm{obs}$" if labeling=="likelihood" else r"$\psi'(s)$" )


    ax.plot ( [C95]*2,      yrange, color=C95_color, linewidth=linewidths, label=r"$\ell_{95}$" if labeling=="likelihood" else r"$C_{95}$" )
    if labeling=="likelihood":
        ax.plot ( [Caic]*2, yrange, color=Caic_color,linewidth=linewidths, label=r"$\ell_\mathrm{AIC}$" )
    elif not Caic is None:
        ax.plot ( [Caic]*2, yrange, color=Caic_color,linewidth=2 )


    if not labeling is "other":
        label_axes ( title="Permutation test", xlabel="log-likelihood" if labeling=="likelihood" else "slope", legend="best", ax=ax )
    ax.set_ylim ( *yrange )

    xti.append ( np.round(observed_value) )
    ax.set_xticks ( xti )

def history_kernels ( estimated_stimulus_kernel, estimated_response_kernel, ci_kernels, ax=None, presentation="left/right", ground_truth=None ):
    """plot history kernels

    :Parameters:
        *estimated_stimulus_kernel*
            stimulus kernel estimated from the data
        *estimated_response_kernel*
            response kernel estimated from the data
        *ci_kernels*
            a sequence of confidence regions for the kernels as returned by
            statistics.history_kernel_ci()
        *ax*
            pylab.axes where the plot should go
        *presentation*
            how should the kernels be presented? Selection of either 'left/right'
            or 'correct/incorrect'

    :Example:
    >>> skernel = [1.2,.5,.3,.1]
    >>> rkernel = [.1,.1,0,0]
    >>> ci_kernels = [ [[1.3,.6,.4,.2],[.8,.3,.1,-.05]],[[.2,.2,.1,.1],[-.05,0.,-.1,-.1]],[[1.5,.8,.5,.3],[.7,.3,0.,-.2]],[[1.2,.5,.5,.2],[.9,.2,0.,-.05]] ]
    >>> history_kernels ( skernel, rkernel, ci_kernels )
    >>> pl.savefig ( 'test/history_kernels.png' ); pl.close()
    """
    if presentation=="left/right":
        kernels = (estimated_stimulus_kernel,estimated_response_kernel)
        colors  = (stimulus_color,response_color)
        labels  = ("stimulus","response")
        if not ci_kernels is None:
            CI      = np.array(ci_kernels[:2])
        else:
            CI      = None
        if not ground_truth is None:
            true_kernels = ground_truth['stimulus_kernel'],\
                    ground_truth['response_kernel']
    elif presentation=="correct/incorrect":
        kernels = (estimated_stimulus_kernel+estimated_response_kernel,-estimated_stimulus_kernel+estimated_response_kernel)
        colors  = (correct_color,incorrect_color)
        labels  = ("correct","incorrect")
        if not ci_kernels is None:
            CI      = np.array(ci_kernels[2:])
        else:
            CI      = None
        if not ground_truth is None:
            true_kernels = ground_truth['stimulus_kernel']+\
                    ground_truth['response_kernel'],\
                    -ground_truth['stimulus_kernel']+\
                    ground_truth['response_kernel']
    else:
        raise ValueError, "presentation should be either 'left/right' or 'correct/incorrect'"

    if CI is None:
        CI = np.array([[kernels[0],kernels[0]],[kernels[1],kernels[1]]])


    if ax is None:
        ax = pl.gca()
    ax = prepare_axes ( ax )

    # Plot confidence regions
    lags = np.arange ( len(estimated_stimulus_kernel) ) + 1
    for i in [0,1]:
        fc = 0.5*np.array(colors[i])+0.5*np.ones(3)
        ax.fill ( np.concatenate ( (lags,lags[::-1]) ), np.concatenate ( (CI[i,0,:],CI[i,1,::-1]) ),
                facecolor=fc, edgecolor=0.5*colors[i], alpha=0.7 )

    kernellines = []
    for i in [0,1]:
        if not ground_truth is None:
            ax.plot ( lags, true_kernels[i], color=0.5*colors[i] )
        kernellines += ax.plot ( lags, kernels[i], 'o',
                    markerfacecolor=colors[i], markeredgecolor=0.5*colors[i], label=labels[i] )

    ax.set_xlim ( 1-0.01*len(estimated_stimulus_kernel),len(estimated_stimulus_kernel)+0.01*len(estimated_stimulus_kernel) )
    ax.set_xticks ( lags )

    # label_axes ( title="history kernels", xlabel="lag", ylabel="equivalent stimulus strength", legend='best', ax=ax )

    return kernellines

def plot_data_summary ( data, ax=None, color='k', label="" ):
    """Plot a data summary with small ellipes

    :Parameters:
        *data*
            an array with three columns: (signed) stimulus intensity, response counts, number of trials
        *ax*
            pylab.axes where the plot should go
        *color*
            color for the datapoints

    :Example:
    >>> import monkey
    >>> d = monkey.MonkeyData ( "../../Data/Steffi/data_from_steffi/dginfokonrad.mat" )
    >>> plot_data_summary ( d.getsummary ( 0 ), color='b', label="Nbw|Sbw" )
    >>> plot_data_summary ( d.getsummary ( 1 ), color='r', label="Nbw|Sc" )
    >>> plot_data_summary ( d.getsummary ( 2 ), color='y', label="Nc|Sbw" )
    >>> pl.savefig ( 'test/montecarlo_test_1.png' ); pl.close()
    >>> import checker
    >>> d = checker.CheckerData ( "../../Data/data_from_marianne/data_sets_fuer_ingo/ip/" )
    >>> plot_data_summary ( d.getsummary ( 0 ), color='b', label="Nbw|Sbw" )
    >>> plot_data_summary ( d.getsummary ( 1 ), color='r', label="Nbw|Sc" )
    >>> plot_data_summary ( d.getsummary ( -1 ), color='y', label="Nc|Sbw" )
    >>> pl.savefig ( 'test/montecarlo_test_2.png' ); pl.close()
    """
    if ax is None:
        ax = pl.gca()
    ax = prepare_axes ( ax )

    x = data[:,0]
    width = 0.008*x.ptp()
    set_label = True

    for s,k,n in data:
        p = float(k)/n
        center = (s,p)
        B = stats.beta ( 1+k,1+n-k )
        a = Arc ( center, width, B.ppf(.975)-p, theta1=0, theta2=180, facecolor=color, edgecolor=color )
        b = Arc ( center, width, B.ppf(.025)-p, theta1=0, theta2=180, facecolor=color, edgecolor=color )
        ax.add_patch ( a )
        ax.add_patch ( b )
        if set_label and not label == "":
            ax.plot ( [s-0.5*width,s+0.5*width], [p]*2, color=color, label=label )
            set_label = False
        else:
            ax.plot ( [s-0.5*width,s+0.5*width], [p]*2, color=color )

    ax.set_ylim ( -.01, 1.01 )

def plot_pmf ( pi, dlal, x=None, ax=None, color='k', alpha=1. ):
    """Plot a psychometric function

    :Example:
    >>> plot_pmf ( [.3,.6,.1], [1.,2.] )
    >>> pl.savefig ( 'test/plot_pmf.png' ); pl.close()
    """
    if ax is None:
        ax = pl.gca ()
    ax = prepare_axes ( ax )

    if x is None:
        x = pl.mgrid[-3:3:100j]

    pmf = ax.plot ( x, pi[1]+pi[2]/(1.+pl.exp(-(dlal[0]+dlal[1]*x))), color=color, linewidth=linewidths, alpha=alpha )
    ax.set_ylim ( -.01, 1.01 )
    return pmf

def plot_nonlinearity_summary ( data, w, pi, ax=None, color='k', label="" ):
    """Plot a data summary after inverting it through the psychometric function
    
    :Example:
    >>> import monkey
    >>> d = pl.array ( zip ( [-2.,-1.,0.,1.,2.],[1,4,10,16,20],[20]*5 ) )
    >>> w = [1.,2.]
    >>> pi = [.07,.05,.85]
    >>> plot_nonlinearity_summary ( d, w, pi )
    >>> pl.savefig ( 'test/plot_nonlinearity_summary.png' )
    """
    def invert ( p ):
        if p<pi[0] or p>1-pi[1]:
            return None
        f = (p-pi[0])/pi[2]
        e = np.log ( f/(1-f) )
        return (e-w[0])/w[1]

    if ax is None:
        ax = pl.gca()
    ax = prepare_axes ( ax )

    x = data[:,0]
    width = 0.008*x.ptp()
    set_label = True

    for s,k,n in data:
        p = invert ( float(k)/n )
        center = (s,p)
        B = stats.beta ( 1+k,1+n-k )
        u = invert ( B.ppf(.975) )
        l = invert ( B.ppf(.025) )
        if p is None:
            continue
        if not u is None:
            a = Arc ( center, width, u-p, theta1=0, theta2=180, facecolor=color, edgecolor=color )
            ax.add_patch ( a )
        if not l is None:
            b = Arc ( center, width, l-p, theta1=0, theta2=180, facecolor=color, edgecolor=color )
            ax.add_patch ( b )
        if set_label and not label == "":
            ax.plot ( [s-0.5*width,s+0.5*width], [p]*2, color=color, label=label )
            set_label = False
        else:
            ax.plot ( [s-0.5*width,s+0.5*width], [p]*2, color=color )

def plot_nonlinearity ( nu, x=None, ax=None, color='k' ):
    """Plot a threshold nonlinearity"""
    from threshold import u_v
    if ax is None:
        ax = pl.gca()
    ax = prepare_axes ( ax )

    if x is None:
        xmin,xmax = ax.get_xlim ()
        x = pl.mgrid[xmin:xmax:100j]

    return ax.plot ( x, u_v ( x, nu ), color=color )
