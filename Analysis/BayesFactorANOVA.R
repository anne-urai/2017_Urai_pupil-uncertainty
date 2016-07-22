# set path
mypath <- '/Users/anne/Data/pupilUncertainty_FigShare'

# activate the BayesFactor package
library("BayesFactor")

# read in data
data = read.csv(sprintf("%s/Data/CSV/ANOVAdat.csv", mypath))
data$subjnr <- factor(data$subjnr)
data$prevPupilBins <- factor(data$prevPupilBins)

# preallocate output
outp = matrix(, 1,2)

# classical rm-ANOVA
an = aov(DV ~ prevPupilBins + Error(subjnr/prevPupilBins), data = data)
sum = summary(an)
outp[1] = sum[[2]][[1]][1,"Pr(>F)"] # WTF is this syntax!

# Bayesian ANOVA
bf = anovaBF(DV ~ prevPupilBins + subjnr, data = data, whichRandom="subjnr")
outp[2] =  exp(bf@bayesFactor$bf) # bf01

# save and get back into mat
colnames(outp) <- c('pvalue', 'bf10')
write.csv(outp, sprintf("%s/Data/CSV/ANOVAresults.csv", mypath))