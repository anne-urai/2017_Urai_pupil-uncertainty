# This code reproduces the analyses in the paper
# Urai AE, Braun A, Donner THD (2016) Pupil-linked arousal is driven 
# by decision uncertainty and alters serial choice bias. 
# 
# Permission is hereby granted, free of charge, to any person obtaining a 
# copy of this software and associated documentation files (the "Software"), 
# to deal in the Software without restriction, including without limitation 
# the rights to use, copy, modify, merge, publish, distribute, sublicense, 
# and/or sell copies of the Software, and to permit persons to whom the 
# Software is furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included 
# in all copies or substantial portions of the Software.
# If you use the Software for your own research, cite the paper.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
# OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
# DEALINGS IN THE SOFTWARE.
#
# Anne Urai, 2016
# anne.urai@gmail.com

# set path
mypath <- '/Users/anne/Data/pupilUncertainty_FigShare'

# activate the BayesFactor package
library("BayesFactor")

# read in data
data = read.csv(sprintf("%s/Data/CSV/sessionANOVAdat.csv", mypath))
data$subjnr <- factor(data$subjnr)
data$prevPupilBins <- factor(data$prevPupilBins)
data$session <- factor(data$session)

# preallocate output

# classical rm-ANOVA
an = aov(DV ~ prevPupilBins*session + Error(subjnr/(prevPupilBins*session)), data = data)
sum = summary(an)
#outp[1] = sum[[2]][[1]][1,"Pr(>F)"] # WTF is this syntax!
#outp[2] = sum[[2]][[1]][1,"F value"] # WTF is this syntax!
#outp[3] = sum[[2]][[1]][1,"Df"] # WTF is this syntax!
#outp[4] = sum[[1]][[1]][1,"Df"] # WTF is this syntax!

outp = matrix(, 1,4)

# Bayesian ANOVA
bf = anovaBF(DV ~ prevPupilBins*session + subjnr, data = data, whichRandom="subjnr")
outp[1:4] = exp(bf@bayesFactor$bf) # bf01

# save and get back into mat
colnames(outp) <- c('pupil', 'session', 'pupil*session', 'subject')
write.csv(outp, sprintf("%s/Data/CSV/sessionANOVAresults.csv", mypath))