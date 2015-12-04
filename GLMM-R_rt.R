######################################################
# load datafile
######################################################

library("lme4")
data <- read.csv("~/Data/UvA_pupil/CSV/2ifc_data2_allsj.csv")

# define equivalent of cirshift in matlab
roll <- function( x , n ){
  if( n == 0 )
    return( x )
  c( tail(x,n) , head(x,-n) )
}

# zscore pupil and motionstrength, log RT
data$decision_pupil <- scale(data$decision_pupil)
data$motionstrength <- scale(data$motionstrength)
data$rt <- scale(log(data$rt)) # log-transform RT

# roll back
data$prevResp <- roll(data$resp, 1) 
data$prevPupil <- roll(data$decision_pupil, 1)
data$prevRT <- roll(data$rt, 1)
data$prevCorr <- roll(data$correct, 1)

# y values must be between 0 and 1
data$resp[data$resp<0] <- 0

# get rid of those trials between blocks
trlcnt <- diff(data$trialnr)
trlcnt <- c(0, trlcnt)
removeTrls <- which(trlcnt != 1)
data <- data[-c(removeTrls), ]

######################################################
# run generalized linear mixed model
######################################################

myMdl.plain = glmer(resp ~ 1 + motionstrength + (1 + motionstrength|subjnr/sessionnr), 
             family = binomial, verbose=0, data = data)
summary(myMdl.plain)
#coefplot2(myMdl.plain)

# now add history terms to the data
myMdl.history = glmer(resp ~ 1 + motionstrength + prevResp + (1 + motionstrength + prevResp|subjnr/sessionnr), 
                   family = binomial, verbose=0, data = data)
summary(myMdl.history)
#coefplot2(myMdl.history)

anova(myMdl.plain, myMdl.history)

######################################################
# include pupil
######################################################

myMdl.pupil = glmer(resp ~ 1 + motionstrength + prevResp + prevPupil + prevResp*prevPupil  
                    + (1 + motionstrength + prevResp + prevPupil + prevResp*prevPupil|subjnr/sessionnr), 
                     family = binomial, verbose=0, data = data)
summary(myMdl.pupil)
anova(myMdl.history, myMdl.pupil) # important, does the pupil model perform better?

######################################################
# separately for correct and error
######################################################

# also make separate regressors for correct and error trials
data$prevRespError <- data$prevResp
data$prevRespError[which(data$prevCorr == 1)] <- 0
data$prevRespCorrect <- data$prevResp
data$prevRespCorrect[which(data$prevCorr == 0)] <- 0

data$prevPupilError <- data$prevPupil
data$prevPupilError[which(data$prevCorr == 1)] <- 0
data$prevPupilCorrect <- data$prevPupil
data$prevPupilCorrect[which(data$prevCorr == 0)] <- 0

# one huge model
myMdl.pupilCorrVsErr = glmer(resp ~ 1 + motionstrength + prevRespCorrect + prevRespError + prevPupilCorrect + prevPupilError 
                    + prevRespCorrect*prevPupilCorrect + prevRespError*prevPupilError
                    + (1 + motionstrength |subjnr/sessionnr), 
                    family = binomial, verbose=0, data = data)
summary(myMdl.pupilCorrVsErr)
#coefplot2(myMdl.pupilCorrVsErr)
anova(myMdl.pupilCorrVsErr, myMdl.pupil)

######################################################
# add additional lags
######################################################
