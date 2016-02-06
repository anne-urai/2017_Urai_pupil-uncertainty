#!/usr/bin/env python

import pylab as pl
import numpy as np
import column # Data set specific
import history,graphics,statistics,model # general
import sys,cPickle
from threshold import u_v
import pdb # debugger
import copy

import glm

__doc__ = """A number of very high level functions


Copyright (C) 2014 Ingo Fruend

This code reproduces the analyses in the paper

    Fruend, Wichmann, Macke (2014): Quantifying the effect of inter-trial dependence on perceptual decisions. J Vis, 14(7): 9.


    Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

    If you use the Software for your own research, cite the paper.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

# Level of logmessages -- 10 allows for INFO but not for DEBUG
#                         20 will suppress INFO
import logging
logging.root.level = 10
logging.BASIC_FORMAT = '%(message)s'

# pylab configuration
pl.rcParams['font.size'] = 9.

####################################### Loading data ##############################

def load_data_file ( filename, header=False, detection=False):
    """Load data set from file and guess initial values
    
    :Parameters:
        *filename*  name of the data file
        *header*    is the first line a header that should be skipped?
        *detection* are the data detection data? For detection data, we
                    fit a threshold nonlinearity on the stimulus values
                    AEU: default = False, since I mainly use discrimination data
        *modulatory* AEU: include the option of having a modulatory term            
    """
    
    h = history.history_impulses ()
    
    # AEU: read text files in a better way
    cdata = np.loadtxt(filename)
    
    # determine whether we have a modulatory factor or not
    nCols = np.shape(cdata)
    if nCols[1] == 6:
            modulatory      = True
    elif nCols[1] == 7: 
        # put two modulatory terms into the same model 
            doublemodulatory = True
    elif nCols[1] == 5:
            modulatory      = False # the default
            doublemodulatory = False
    
    # make the design matrix out of this data
    data = column.ColumnData (
           cdata, impulse_responses=h, threshold=detection, ground_truth=None, modulation=modulatory)

    conditions = np.unique ( cdata[:,1] )

    w0 = np.zeros ( len(conditions)+1, 'd' )
    for c in xrange ( len(conditions) ):
        d = np.array ( data.getsummary ( c ) )
        w = history.get_w0 ( d )
        w0[0] += w[0]
        w0[1+c] = w[1]
    w0[0] /= len(conditions)

    # AEU: changed colors to make them sequential across sessions
    color_list = pl.cm.pink(np.linspace(0, 1, len(conditions)))
    color_list = color_list.tolist() # change to list format
    plotinfo = {
            'colors': color_list,
            'labels': ['1','2','3','4','5','6'],
            'conditions': range(len(conditions)),
            'indices': [[0,1],[0,2],[0,3],[0,4],[0,5],[0,6]],
            'xmin': data.X[:,1].min(),
            'xmax': data.X[:,1].max()
            }

    data.detection = detection
    return data,w0,plotinfo

############################## analyze data ##############################

def search_for_start ( r, X, w0, applythreshold, hf0, pm=(.85,.93,1.,1.07,1.15), storeopt=False, modulation=False):
    """Search for starting values

    :Parameters:
        *d*
            data set
        *w0*
            coarse starting values
        *pm*
            increments or decrements
            
        # AEU ToDo: if d.h, make sure this spits out 3 models, one without history, one with history but without modulation, and one with both
    """

    logging.info ( "Initial preoptimization" )
    Mwh = model.history_model ( r, X[:,:hf0],
            applythreshold=applythreshold,
            w0=w0,
            lm=0.1, hf0=hf0, emiter=100 )
    Mnh = model.history_model ( r, X[:,:hf0],
            applythreshold=applythreshold,
            w0=w0,
            lm=0.1, hf0=hf0, emiter=100 )  
    M0  = model.history_model ( r, X[:,:hf0],
            applythreshold=applythreshold,
            w0=w0,
            lm=0.1, hf0=hf0, emiter=100 )
            
    if modulation:        
        Mhmod = model.history_model ( r, X[:,:hf0],
            applythreshold=applythreshold,
            w0=w0,
            lm=0.1, hf0=hf0, emiter=100 )
            
        # add the history weights the model with history
        Mwh.w = np.concatenate ( (Mwh.w,np.zeros(6,'d')) )
        Mwh.X = copy.copy(X) 
        Mwh.X = Mwh.X[:,:-2] # only the part that has no modulation
        
        # add those same weights + modulatory terms to the hmod model
        Mhmod.w = np.concatenate ( (Mhmod.w,np.zeros(15,'d')) )
        Mhmod.X = X # here, something goes wrong!
        
        # to check that the sizes work
        print "X",np.shape(X)
        print "Mwh",Mwh.w,np.shape(Mwh.X)
        print "Mhmod",Mhmod.w,np.shape(Mhmod.X)
        
    else:
        Mwh.w = np.concatenate ( (Mwh.w,np.zeros(6,'d')) )
        Mwh.X = X 
        Mhmod = [] # return empty

    nhind = 0
    whind = 0
    i = 1
    for al in pm:
        for lm in pm:
            logging.info ( "::::: Optimizing from starting value %d :::::" % (i,) )
            w0 = M0.w.copy()
            w0[1:hf0] *= al
            p0 = M0.pi.copy()
            p0[0] *= lm;
            p0[-1] = 1-p0[0]-p0[1]

            if modulation:
                M_ = model.history_model ( r, X,
                        applythreshold=applythreshold,
                        w0=w0, p0=p0, nu0=M0.nu,
                        lm=0.1, hf0=hf0, verbose=True, emiter=300, storeopt=storeopt )
                if Mhmod.loglikelihood < M_.loglikelihood:
                    logging.info ( "  *model chosen for history + modulation*" )
                    Mhmod = M_
                    whind = i

                M_ = model.history_model ( r, X[:,:-2],
                    applythreshold=applythreshold,
                    w0=w0, p0=p0, nu0=M0.nu,
                    lm=0.1, hf0=hf0, verbose=True, emiter=300, storeopt=storeopt )
                if Mwh.loglikelihood < M_.loglikelihood:
                    logging.info ( "  *model chosen for history*" )
                    Mwh = M_
                    whind = i
            else:
                M_ = model.history_model ( r, X,
                    applythreshold=applythreshold,
                    w0=w0, p0=p0, nu0=M0.nu,
                    lm=0.1, hf0=hf0, verbose=True, emiter=300, storeopt=storeopt )
                if Mwh.loglikelihood < M_.loglikelihood:
                    logging.info ( "  *model chosen for history*" )
                    Mwh = M_
                    whind = i        
                
            M_ = model.history_model ( r, X[:,:hf0],
                    applythreshold=applythreshold,
                    w0=w0, p0=p0, nu0=M0.nu,
                    lm=0.1, hf0=hf0, verbose=True, emiter=300, storeopt=storeopt )
            if Mnh.loglikelihood < M_.loglikelihood:
                logging.info ( "  *model chosen for independent*" )
                Mnh = M_
                nhind = i
            i += 1
            
    logging.info ( "Mwh.w = %s\nMnh.w = %s" % (str(Mwh.w),str(Mnh.w)) )
    logging.info ( "Mwh.ll = %g\nMnh.ll = %s" % (Mwh.loglikelihood,Mnh.loglikelihood) )
    logging.info ( "Starting values:\n  with history: %d\n  without history: %d\n" % (whind,nhind) )

    # NOW, THE HISTORY ONLY MODEL HAS SIZE 22!  
    print "X",np.shape(X)
    print "Mwh",Mwh.w,np.shape(Mwh.X)
    if modulation:
        print "Mhmod",Mhmod.w,np.shape(Mhmod.X)    

    return Mnh,Mwh,Mhmod

def analysis ( d, w0, nsamples=200, perm_collector=statistics.EvaluationCollector ):
    """Analyze a dataset

    :Parameters:
        *d*
            a history.DataSet instance (typically a subclass of history.DataSet)
        *w0*
            starting values for the first parameters in the model. The remaining parameters will
            start at 0. It seems to be a good idea to give starting values for the stimulus dependent
            parameters in the model and leave the parameters for history features at 0.
        *nsamples*
            number of samples for monte carlo procedures
        *perm_collector*
            the collector object for the permutation tests. The default one should do for all
            experiments with a design matrix in which the first column is 1 and the second column
            refers to the slope.

    :Example:
    >>> d,w,plotinfo = load_plaid ()
    >>> results = analysis ( d, w, 10 )
    >>> results.keys()
    ['model_nohist', 'model_w_hist', 'bootstrap', 'permutation_nh', 'permutation_wh']
    """
    
    # get a version of the data without history and without modulation
    if d.__dict__.has_key ( 'detection' ) and d.detection:
        dnh = d.__class__ ( d.fname, threshold=len(d.th_features)>0 )
    else:
        dnh = d.__class__( d.fname, ) # reload the data without history
        dnh.modulation = False
        
    if d.modulation:
        dhmod = copy.copy(d) # this includes everything
        dhmod.modulation = True
        
        # the data with history, but without modulation
        dwh = copy.copy(d)
        dwh = dwh.__class__( dwh.fname, dwh.h, False, None, False ) 
        dwh.modulation = False
    else:
        dwh = d
      
    # return indices of difficult trials (p > 0.75) and difficult trials (p < 0.55)
    easy,difficult = d.performance_filter ()
    logging.info ( "Fitting models" )
        
    # check if we have the right sizes here
    print "d", np.shape(d.X), np.shape(d.r)
    print "dnh", np.shape(dnh.X), np.shape(dnh.r)
    print "dwh", np.shape(dwh.X), np.shape(dwh.r)
    if d.modulation:
        print "dhmod", np.shape(dhmod.X), np.shape(dhmod.r)
    # here, the sizes seem OK

    Mnh,Mwh,Mhmod = search_for_start ( d.r, d.X, w0, d.th_features, d.hf0, storeopt=True, modulation=d.modulation)
    logging.info ( "likelihood for independent responses: %g" % (Mnh.loglikelihood,) )
    logging.info ( "likelihood for history model: %g" % (Mwh.loglikelihood,) )
    if d.modulation:    
        logging.info ( "likelihood for modulation + history model: %g" % (Mhmod.loglikelihood,) )

    # here, the Mwh.X has the wrong shape!
    print "nh",Mnh.w,Mnh.pi,np.shape(Mnh.X)
    print "wh",Mwh.w,Mwh.pi,np.shape(Mwh.X)
    if d.modulation:
        print "hmod",Mhmod.w,Mhmod.pi,np.shape(Mhmod.X)

    print 'start monte carlo'
    # Monte Carlo testing
    if nsamples>0:
        
        # now, dhmod has too few columns
        
        if d.modulation:
            print 'permuting with modulation'
            r_,X_ = dhmod.permutation () # permute the whole thing
        else:
            r_,X_ = dwh.permutation () # permute the whole thing

        Mnh_perm, Mwh_perm, Mhmod_perm = search_for_start ( r_, X_,
                w0, d.th_features, d.hf0, modulation=d.modulation)

        # Set the states of the random number generators to the same values to get the exact same sequence of random numbers
        dnh.rng.set_state ( dwh.rng.get_state () ) 
        if d.modulation:
            dhmod.rng.set_state ( dwh.rng.get_state () ) 
    
        logging.info ( "Permutation without history" )
        perm_collector = statistics.EvaluationCollector (
               Mnh, easy, difficult )
        permutation_nh   = pl.array ( statistics.mcsimulation (
            dnh.permutation, perm_collector,
            nsamples, Mnh_perm.w, Mnh_perm.pi, Mnh_perm.nu,
            verbose=logging.root.level<20,
            hf0=dnh.hf0, applythreshold=Mwh.applythreshold ) )
        
        logging.info ( "Permutation with history, no modulation" )
        perm_collector = statistics.EvaluationCollector (
                Mwh, easy, difficult )
        permutation_wh = pl.array ( statistics.mcsimulation (
            dwh.permutation, perm_collector,
            nsamples, Mwh_perm.w, Mwh_perm.pi,
            Mwh_perm.nu, verbose=logging.root.level<20,
            hf0=dwh.hf0, applythreshold=Mwh.applythreshold ) )
       
        if d.modulation:
            logging.info ( "Permutation with history and modulatory interaction" )
            perm_collector = statistics.EvaluationCollector (
                Mhmod, easy, difficult )

                # this is where the bug occurs!
                # do we have enough info about dhmod.permutation?
            print(dhmod.modulation)
            permutation_hmod = pl.array ( statistics.mcsimulation (
                dhmod.permutation, perm_collector,
                nsamples, Mhmod_perm.w, Mhmod_perm.pi,
                Mhmod_perm.nu, verbose=logging.root.level<20,
                hf0=dhmod.hf0, applythreshold=Mhmod.applythreshold ) )

        logging.info ( "Bootstrap" )
        if d.modulation:
            kcollector = statistics.Kernel_and_Slope_Collector (
                    dhmod.h, dhmod.hf0, slopeindex=range(1,dhmod.hf0) )
            bootstrap   = pl.array ( statistics.mcsimulation (
                dhmod.bootstrap, kcollector,
                nsamples, Mhmod.w, Mhmod.pi, Mhmod.nu,
                verbose=logging.root.level<20,
                hf0=dhmod.hf0, applythreshold=Mhmod.applythreshold ) )
        else:
            # only bootstrap without the modulatory interaction
            kcollector = statistics.Kernel_and_Slope_Collector (
                    dwh.h, dwh.hf0, slopeindex=range(1,dwh.hf0) )
            bootstrap   = pl.array ( statistics.mcsimulation (
                dwh.bootstrap, kcollector,
                nsamples, Mwh.w, Mwh.pi, Mwh.nu,
                verbose=logging.root.level<20,
                hf0=dwh.hf0, applythreshold=Mwh.applythreshold ) )
                
    else:
        permutation_wh  = None
        permutation_nh  = None
        permutation_hmod = None
        bootstrap       = None

    if d.modulation:
        results = {
            'model_nohist': Mnh,
            'model_w_hist': Mwh,
            'model_h_mod': Mhmod,
            'permutation_wh': permutation_wh,
            'permutation_nh': permutation_nh,
            'permutation_hmod': permutation_hmod,
            'bootstrap': bootstrap
            }
    else:
        results = {
            'model_nohist': Mnh,
            'model_w_hist': Mwh,
            'permutation_wh': permutation_wh,
            'permutation_nh': permutation_nh,
            'bootstrap': bootstrap
            }

    return results

############################## Display data ##############################

def plot ( d, results, infodict ):
    """plot all results"""

    fig = pl.figure ()
    ax  = graphics.canonical_axes ( fig )

    required_infokeys = ['labels','colors','indices','conditions','xmin','xmax']
    for k in required_infokeys:
        if not k in infodict.keys():
            raise ValueError, "Key %s was not in infodict" % (k,)

    pmfplot ( d, results, infodict, ax.pmf )
    nonlinearityplot ( d, results, infodict, ax.uv )
    permutationplot ( d, results, infodict, ax.likeli )
    kernelplot ( d, results, infodict, ax.history_rz, ax.history_perf )
    slopeplot ( d, results, infodict, ax.slopes )
    
    # AEU: add an extra row of plots with the modulatory pupil weights
    

############################## high level plotting routines -- called internally by plot()

def pmfplot ( d, results, infodict, ax, errors=True ):
    """Generate the pmf plot"""
    for i,c in enumerate(infodict['conditions']):
        c = int(c)
        d_ = d.getsummary ( c )
        x = pl.mgrid[infodict['xmin']:infodict['xmax']:100j]
        if len(d.th_features)>0:
            # d_[:,0] = u_v ( d_[:,0], results['model_w_hist'].nu )
            # x = u_v ( x, results['model_w_hist'].nu )
            d_[:,0] = u_v ( d_[:,0], results['model_nohist'].nu )
            x = u_v ( x, results['model_nohist'].nu )

        if errors:
            graphics.plot_data_summary ( d_, ax,
                infodict['colors'][i], infodict['labels'][i] )
        else:
            # AEU: change the marker to make datapoints small. What goes wrong in the plot here?
            ax.plot ( d_[:,0], d_[:,1] / d_[:,2], '.', color=infodict['colors'][i], label=infodict['labels'][i] )

        # wfit  = results['model_w_hist'].w[infodict['indices'][c]]
        wfit  = results['model_nohist'].w[infodict['indices'][i]]
        w0fit = results['model_w_hist'].w0[infodict['indices'][i]]
        pfit  = results['model_w_hist'].pi
        p0fit  = results['model_nohist'].pi

        if not d.ground_truth is None:
            wgfit = d.ground_truth['w'][infodict['indices'][i]]
            pgfit = d.ground_truth['pi']
            gt = graphics.plot_pmf ( pgfit, wgfit, x, ax,
                    (np.array([255,240,240],'d')/255+infodict['colors'][i])/2. )
            pl.setp (gt, linestyle='--' )
 
        # graphics.plot_pmf ( pfit, w0fit, x, ax, [.9,.9,.9], alpha=0.1 )
        graphics.plot_pmf ( p0fit, wfit,  x, ax, infodict['colors'][i] )
    graphics.label_axes (
            title="(A) psychometric function",
            xlabel=r"transduced stimulus $u_\nu(s\tilde{z})$",
            ylabel=r"probability for $r=1$",
            nxticks=5,
            ax=ax )

def nonlinearityplot ( d, results, infodict, ax ):
    """Plot with the nonlinearity"""
    xmin,xmax = infodict['xmin'],infodict['xmax']
    M = results['model_w_hist']
    for c in infodict['conditions']:
        d_ = d.getsummary ( c )

        graphics.plot_nonlinearity_summary ( d_, M.w, M.pi, ax, color=infodict['colors'][c], label=infodict['labels'][c] )
    graphics.plot_nonlinearity ( M.nu, np.mgrid[xmin:xmax:100j], ax, 'k' )
    if not d.ground_truth is None:
        nl = graphics.plot_nonlinearity ( d.ground_truth['nu'], np.mgrid[xmin:xmax:100j], ax, np.array([255,230,230],'d')/255 )
        pl.setp ( nl, linestyle='--' )
    graphics.label_axes ( title="(B) nonlinearity",
            xlabel=r"raw stimulus $s\tilde{z}$",
            ylabel=r"transduced stimulus $u_\nu(s\tilde{z})$",
            nxticks=5,
            ax=ax)

def permutationplot ( d, results, infodict, ax, noaic=False ):
    """permutation test"""
    l_obs = results['model_w_hist'].loglikelihood
    Caic  = results['model_nohist'].loglikelihood+len(results['model_w_hist'].w)-len(results['model_nohist'].w)
    out = [l_obs]
    print "l_obs=",l_obs
    print "Caic= ",Caic
    if not results['permutation_wh'] is None:
        hist, C95 = statistics.historytest ( results['permutation_wh'][:,0] )
        print "C95=  ",C95
        out.append ( C95 )
        out.append ( Caic )
        out.append ( np.mean ( results['permutation_wh'][:,0]<l_obs ) )
    else:
        hist = None
        C95  = None
        out.append ( None )
        out.append ( Caic )
        out.append ( None )
    if noaic:
        Caic=None
    graphics.montecarlo_test ( l_obs, hist, C95, Caic, ax, "likelihood" )
    graphics.label_axes ( title="(C) Permutation test",
            xlabel="log-likelihood",
            ax=ax )
    return out

def kernelplot ( d, results, infodict, ax1, ax2, legend='lower right' ):
    """Plot historykernels"""
    M = results['model_w_hist']
    bootstrap = results['bootstrap']

    C = statistics.Kernel_and_Slope_Collector ( d.h, d.hf0, range(1,d.hf0) )
    K = C(M)
    print 'K', K
    print d.hf0
    nlags = d.h.shape[0]
    print 'nlags = ', nlags
   
    hr = K[:nlags]
    hz = K[nlags:2*nlags]
    # what is stored here?
    hr *= K[-2]
    hz *= K[-2]   
    
    print "h_r[1]",hr[0]
    print "h_z[1]",hz[0]

    if bootstrap is None:
        kernellen = (bootstrap.shape[1]-2)/2
        print kernellen
        al = bootstrap[:,-2]
        al.shape = (-1,1)
        bootstrap[:,:-2] *= al # Like that?
        print K[-2],pl.prctile(bootstrap[:,-2]),pl.mean(bootstrap[:,-2])

        hci = statistics.history_kernel_ci (
                bootstrap[:,kernellen:-2], bootstrap[:,:kernellen],
                hz, hr )
    else:
        hci = None

    kl = graphics.history_kernels ( hz, hr, hci, ax1,   "left/right", ground_truth = d.ground_truth )
    kl += graphics.history_kernels ( hz, hr, hci, ax2, "correct/incorrect", ground_truth = d.ground_truth )

    labely,labelh = same_y ( ax1, ax2 )

    graphics.label_axes ( title="(D) stimulus and response kernels",
            xlabel="lag",
            ylabel="equivalent stimulus strength",
            legend=legend,
            ax=ax1 )
    graphics.label_axes ( title="(E) correct and incorrect kernels",
            xlabel="lag",
            ylabel="equivalent stimulus strength",
            legend=legend,
            ax=ax2 )
    # pl.setp ( (ax1,ax2), ylim=(-6,6) )

    return kl
    
def kernelplot_Mod ( d, results, infodict, ax1, ax2, ax3, legend='lower right' ):
    """Plot historykernels"""
    M = results['model_w_hist']
    bootstrap = results['bootstrap']

    C = statistics.Kernel_and_Slope_Collector ( d.h, d.hf0, range(1,d.hf0) )
    K = C(M)
    print 'K', K
    print d.hf0
    print 'd.h.shape[0]', d.h.shape[0]
    
    print 'adding modulation kernels'
    hr = K[:d.h.shape[0]]
    hz = K[d.h.shape[0]:2*d.h.shape[0]]
    hr_pupil = K[2*d.h.shape[0]:3*d.h.shape[0]]
    hz_pupil = K[3*d.h.shape[0]:-2] 
    hr_pupil *= K[-2]
    hz_pupil *= K[-2]
   
    if bootstrap is None:
        kernellen = (bootstrap.shape[1]-2)/2
        print kernellen
        al = bootstrap[:,-2]
        al.shape = (-1,1)
        bootstrap[:,:-2] *= al # Like that?
        print K[-2],pl.prctile(bootstrap[:,-2]),pl.mean(bootstrap[:,-2])

        hci = statistics.history_kernel_ci (
                bootstrap[:,kernellen:-2], bootstrap[:,:kernellen],
                hz, hr )
    else:
        hci = None

    kl = graphics.history_kernels ( hz, hr, hci, ax1,   "left/right", ground_truth = d.ground_truth )
    kl += graphics.history_kernels ( hz, hr, hci, ax2, "correct/incorrect", ground_truth = d.ground_truth )

    labely,labelh = same_y ( ax1, ax2 )

    graphics.label_axes ( title="(D) stimulus and response kernels",
            xlabel="lag",
            ylabel="equivalent stimulus strength",
            legend=legend,
            ax=ax1 )
    graphics.label_axes ( title="(E) correct and incorrect kernels",
            xlabel="lag",
            ylabel="equivalent stimulus strength",
            legend=legend,
            ax=ax2 )
    # pl.setp ( (ax1,ax2), ylim=(-6,6) )

    return kl

def slopeplot ( d, results, infodict, ax ):
    """slope results of the permutation test"""

    ax = graphics.prepare_axes ( ax, haveon=('bottom',) )

    h = np.histogram ( results['permutation_wh'][:,1] )
    graphics.montecarlo_test ( results['model_w_hist'].w[1], h, pl.prctile(results['permutation_wh'][:,1], 95), ax=ax, labeling='slope' )
    # ax.bar ( b[:-1], h, pl.diff(b), edgecolor=graphics.histogram_color, facecolor=graphics.histogram_color )
    # yrange = ax.get_ylim ()
    # ax.axvline ( results['model_w_hist'].w[1], ymin=yrange[0], ymax=yrange[1] )
    # ax.set_ylim ( yrange )
    # # ax.set_xlim ( .3,3 )
    # # ax.set_xticks ( [.5,1,1.5,2,2.5,3] )

    graphics.label_axes ( title="(F) slope effects",
            xlabel=r"slope",
            ax = ax )

def same_y ( *axes ):
    """Set all axes to the same ylimits"""
    lims = []
    for ax in axes:
        lims.append ( ax.get_ylim () )
    lims = pl.array(lims)
    ymin = abs(lims.min())
    ymax = abs(lims.max())
    yl = max(ymin,ymax)
    for ax in axes:
        ax.set_ylim ( -yl, yl )
    return -.8*yl,.1*yl
