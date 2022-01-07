#!/usr/bin/env python

__doc__ = """

Copyright (C) 2014 Ingo Fruend

This code reproduces the analyses in the paper

    Fruend, Wichmann, Macke (2014): Quantifying the effect of inter-trial dependence on perceptual decisions. J Vis, 14(7): 9.


    Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

    If you use the Software for your own research, cite the paper.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

import history
import numpy as np

class ColumnData ( history.DataSet ):
    def __init__ ( self, data, impulse_responses=None, threshold=False, ground_truth=None, modulation=False, doublemodulation=False):
        """A data set consisting of multiple columns of data

        :Parameters:
            *data*
                an array with 5 columns (block,condition,stimulus,target,response)
                block should be positive integers, condition, should be positive
                integers, stimulus should be positive, target should have values of
                0 and 1, response should have values of 0 and 1.
                AEU: data can contain a 6th column, indicating a modulatory term
            *impulse_responses*
                an array with the impulse responses of the history filters in the
                columns. Such an array is most easily constructed using the function
                history.history_impulses.
            *threshold*
                set this to True, if you want the stimulus to be thresholded.
            *ground_truth*
                for simulated data, this can be the model instance that
                contains the generating parameters

        :Example:
        >>> c = np.array ( [[1,1, 1,0,1], [1,1, 1,1,1], [1,1, 1,0,0], [2,1,.5,1,1], [2,1,.5,0,1], [3,1,.3,1,0]] )

        Example without thresholding
        >>> d = ColumnData ( c, None )
        >>> d.X
        array([[ 1. , -1. ],
               [ 1. ,  1. ],
               [ 1. , -1. ],
               [ 1. ,  0.5],
               [ 1. , -0.5],
               [ 1. ,  0.3]])
        >>> d.r
        array([ 1.,  1.,  0.,  1.,  1.,  0.])
        >>> d.th_features
        []
        >>> d.hf0
        2
        >>> d.getsummary ()
        array([[-1. ,  1. ,  2. ],
               [ 1. ,  1. ,  1. ],
               [-0.5,  1. ,  1. ],
               [ 0.5,  1. ,  1. ],
               [ 0.3,  0. ,  1. ]])


        Example with thresholding
        >>> d_th = ColumnData ( c, [], True )
        >>> d_th.th_features
        [1]

        Example with multiple conditions
        >>> c = np.array ( [[1,1, 1,0,1],[1,1, 1,1,1],[2,2, 1,1,0],[2,2, 1,0,1],[3,1,.5,0,0],[3,1,.5,1,0]] )
        >>> d_m = ColumnData ( c, [] )
        >>> d_m.X
        array([[ 1. , -1. ,  0. ],
               [ 1. ,  1. ,  0. ],
               [ 1. ,  0. ,  1. ],
               [ 1. ,  0. , -1. ],
               [ 1. , -0.5,  0. ],
               [ 1. ,  0.5,  0. ]])
        >>> d_m.hf0
        3
        """
        history.DataSet.__init__ ( self, impulse_responses )
        self.__blocks     = np.sort ( np.unique ( data[:,0] ) )
        self.__conditions = np.sort ( np.unique ( data[:,1] ) )
        self.__data       = data
        self.fname        = data
        self.modulation   = modulation
        self.doublemodulation   = doublemodulation
        self.__construct_design ( )
        self.__threshold  = threshold
        self.__th_features = []
        self.ground_truth = ground_truth

        if threshold:
            for i,condition in enumerate ( self.__conditions ):
                self.__th_features.append ( 1+i )
    @property
    def X ( self ):
        """Design matrix"""
        return self.__X
    @property
    def r ( self ):
        """response vector"""
        return self.__r

    def permutation ( self ):

        if self.modulation:
            """Permute only the pupil part"""
            data = self.__data.copy()
            for condition in self.__conditions:
                cond_idx = data[:,1] == condition
                these_data = data[cond_idx,-1] # get only the pupil values
                np.random.shuffle(these_data)
                data[cond_idx,-1] = these_data

            # then construct designM again
            C = ColumnData ( data, self.h, self.__threshold, ground_truth=None, modulation=self.modulation )
            return C.r,C.X

        elif self.doublemodulation:
            """Permute both of the modulatory regressors independently of one another"""
            data = self.__data.copy()
            for condition in self.__conditions:
                cond_idx = data[:,1] == condition
                these_data = data[cond_idx,-2:] # get the pupil and RT values
                np.random.shuffle(these_data)
                data[cond_idx,-2:] = these_data

            # then construct designM again
            C = ColumnData ( data, self.h, self.__threshold, ground_truth=None, doublemodulation=self.doublemodulation )
            return C.r,C.X

        else:
            """Return a conditionwise permutation of the original dataset"""

            data = self.__data.copy()
            for condition in self.__conditions:
                cond_idx = data[:,1] == condition
                these_data = data[cond_idx,1:] # get the data within a session, leave blocknrs intact
                np.random.shuffle(these_data)
                data[cond_idx,1:] = these_data

            # for block in self.__blocks:
            #    block_index = self.__data[:,0] == block
            #    these_data  = self.__data[block_index,:]
            #    np.random.shuffle ( these_data )
            #    data[block_index,:] = these_data

            # after shuffling all trials within a block, recompute the designM
            # print((self.h))
            C = ColumnData ( data, self.h, self.__threshold, ground_truth=None, modulation=self.modulation)

            return C.r,C.X

    @property
    def th_features ( self ):
        """Features to be thresholded"""
        return self.__th_features

    @property
    def hf0 ( self ):
        """Starting index of history features"""
        return 1+len(self.__conditions)

    def getsummary ( self, condition=0 ):
        """A three column summary of the data from one condition"""
        condition = self.__conditions[condition]
        condition_index = self.__data[:,1]==condition
        these_data = self.__data[condition_index,:]
        blocks = np.unique ( these_data[:,0] )
        out = []
        for block in blocks:
            block_index = these_data[:,0] == block
            block_data = these_data[block_index,:]
            stimuli    = block_data[:,2]*np.array( [history.get_code ( z, [-1,1],[0,1] ) for z in block_data[:,3]] )
            for stimulus in np.unique (stimuli):
                stim_index = stimuli == stimulus
                r = block_data[stim_index,4].sum()
                n = len(block_data[stim_index,4])
                out.append ( (stimulus,r,n) )
        return np.array ( out )

    def __construct_design ( self ):
        """Construct the design matrix
        # AEU: add another column for a modulatory weight
        # contingent on having the -modulatory option in self
        """
        x = []
        y = []
        p = []
        codes_z = [0,1]
        codes_r = [0,1]
        nconditions = len(self.__conditions)

        for block in self.__blocks:
            block_index = self.__data[:,0] == block
            ntrials_this_block = block_index.sum()
            these_data         = self.__data[block_index,:]

            z = these_data[:,3]
            r = these_data[:,4]
            d = these_data[:,2] # coherence levels

            if self.modulation:
                pupil = these_data[:,5] # modulatory
            elif self.doublemodulation:
                pupil = these_data[:,5] # modulatory 1, in this case pupil
                reactiontime = these_data[:,6] # modulatory 2, in this case rt

            x_ = np.zeros ( (ntrials_this_block,1+nconditions) )
            x_[:,0] = 1.

            for i,condition in enumerate ( self.__conditions ):
                condition_index = these_data[:,1] == condition
                condition_data  = these_data[condition_index,:]
                x_[condition_index,1+i] = condition_data[:,2]*(2*condition_data[:,3]-1)

            z_ = np.zeros ( z.shape )
            r_ = np.zeros ( z.shape )

            if self.modulation: # initialize main and interaction terms as well
                p_ = np.zeros ( z.shape )
                p_z = np.zeros ( z.shape )
                p_r = np.zeros ( z.shape )
            elif self.doublemodulation: # initialize main and interaction terms for both modulators
                p_ = np.zeros ( z.shape )
                p_z = np.zeros ( z.shape )
                p_r = np.zeros ( z.shape )
                rt_ = np.zeros ( z.shape )
                rt_z = np.zeros ( z.shape )
                rt_r = np.zeros ( z.shape )

            # loop over the nr of lags
            for i in range ( len(z) ):
                z_[i] = history.get_code ( z[i], [-1,1], codes_z )
                r_[i] = history.get_code ( r[i], [-1,1], codes_r )

                if self.modulation:
                    p_[i] = pupil[i] # no coding needed, this is a continuous measure
                    p_z[i] = z_[i]*pupil[i] # interaction per lag
                    p_r[i] = r_[i]*pupil[i]
                elif self.doublemodulation:
                    p_[i] = pupil[i] # no coding needed, this is a continuous measure
                    p_z[i] = z_[i]*pupil[i] # interaction per lag
                    p_r[i] = r_[i]*pupil[i]
                    rt_[i] = reactiontime[i] # no coding needed, this is a continuous measure
                    rt_z[i] = z_[i]*reactiontime[i] # interaction per lag
                    rt_r[i] = r_[i]*reactiontime[i]

            # convolve each set of nlag regressors with the exponential history kernels
            hr = history.history_features ( self.h, r_ )
            hz = history.history_features_stim ( self.h, z_, d )

            if self.modulation:
                # main and interaction terms as well
                hp = history.history_features ( self.h, p_ )
                hpr = history.history_features ( self.h, p_r )
                hpz = history.history_features_stim ( self.h, p_z, d )
            elif self.doublemodulation:
                # main and interaction terms as well
                hp = history.history_features ( self.h, p_ )
                hpr = history.history_features ( self.h, p_r )
                hpz = history.history_features_stim ( self.h, p_z, d )
                hrt = history.history_features ( self.h, rt_ )
                hrtr = history.history_features ( self.h, rt_r )
                hrtz = history.history_features_stim ( self.h, rt_z, d )

            # now concatenate all columns into a design matrix
            if not hr is None:
                x_ = np.c_[x_,hr]
            if not hz is None:
                x_ = np.c_[x_,hz]

            # append all the interaction regressors
            if self.modulation:
                if not hp is None:
                    x_ = np.c_[x_,hp]
                if not hpr is None:
                    x_ = np.c_[x_,hpr]
                if not hpz is None:
                    x_ = np.c_[x_,hpz]
            if self.doublemodulation:
                if not hp is None:
                    x_ = np.c_[x_,hp]
                if not hpr is None:
                    x_ = np.c_[x_,hpr]
                if not hpz is None:
                    x_ = np.c_[x_,hpz]
                if not hrt is None: # also the RT regressors
                    x_ = np.c_[x_,hrt]
                if not hrtr is None:
                    x_ = np.c_[x_,hrtr]
                if not hrtz is None:
                    x_ = np.c_[x_,hrtz]

            # how often was the subject correct?
            correct = z_==r_
            performance = np.mean ( correct )

            y.append ( r )
            x.append ( x_ )
            p.append ( performance+np.zeros ( r.shape ) )

        self.__X = np.concatenate ( x, 0 )
        self.__r = np.concatenate ( y, 0 )
        self.__p = np.concatenate ( p, 0 )

    def performance_filter ( self, p0=0.75, p1=0.55 ):
        """Return indices of easy trials and difficult trials

        easy trials: performance better than p1
        difficult trials: performance between p1 and p0
        """

        easy = self.__p>p0
        difficult = np.logical_and ( self.__p>p1, self.__p<=p0 )
        return easy, difficult

    @property
    def design ( self ):
        """the design used by the constructor"""
        return self.__data.copy()

if __name__ == "__main__":
    import doctest
    doctest.testmod ()
