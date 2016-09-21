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