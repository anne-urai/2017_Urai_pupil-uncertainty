# set path
mypath <- '/Users/anne/Data/pupilUncertainty_FigShare'

# activate the BayesFactor package
library("BayesFactor")

# read in data
data = read.csv(sprintf("%s/Data/CSV/ANOVAdat.csv", mypath))
data$subjnr <- factor(data$subjnr)
data$prevPupilBins <- factor(data$prevPupilBins)

# preallocate output
outp = matrix(, 1,5)

# classical rm-ANOVA
an = aov(DV ~ prevPupilBins + Error(subjnr/prevPupilBins), data = data)
sum = summary(an)
outp[1] = sum[[2]][[1]][1,"Pr(>F)"] # WTF is this syntax!
outp[2] = sum[[2]][[1]][1,"F value"] # WTF is this syntax!
outp[3] = sum[[2]][[1]][1,"Df"] # WTF is this syntax!
outp[4] = sum[[1]][[1]][1,"Df"] # WTF is this syntax!

# Bayesian ANOVA
bf = anovaBF(DV ~ prevPupilBins + subjnr, data = data, whichRandom="subjnr")
outp[5] =  exp(bf@bayesFactor$bf) # bf01

# save and get back into mat
colnames(outp) <- c('pvalue', 'F', 'df1', 'df2', 'bf10')
write.csv(outp, sprintf("%s/Data/CSV/ANOVAresults.csv", mypath))