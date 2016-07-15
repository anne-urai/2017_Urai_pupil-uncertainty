# set path
mypath <- '/Users/anne/Data/pupilUncertainty_FigShare'


# ============================================ #
# LAVAAN PACKAGE
# ============================================ #

# activate the lavaan package
library("lavaan", lib.loc="~/Library/R/3.2/library")

# preallocate
nr_subjects = 27;
data = matrix(, nrow = nr_subjects, ncol = 7)

for ( s in 1:27 ) {
  
  # load data that have been coded to indicate repetition
  mydata = read.csv(sprintf("%s/Data/CSV/SEMdata_sj%02d.csv", mypath, s))
  
  # make sure repeat is a logical array, so logistic regression is used
  mydata$repeat. <- factor(mydata$repeat.)
  
  # specify the naive model
  model <- ' # direct effect
			      repeat. ~ c*uncertainty + cov1*stimrepeat
          '
  fit <- sem(model, data=mydata, ordered=c("repeat."))
  coefs <- parameterEstimates(fit)
  data[s, 1] = coefs[['est']][1]
  
  # specify the mediation model
  model <- '
            # indirect effect through mediator
            decision_pupil ~ a*uncertainty
            repeat. ~ b*decision_pupil + c1*uncertainty + cov2*stimrepeat

            # specify the effects we want to test
            indirect := a*b
            direct   := c1
            total    := c1 + (a*b)
            '
  
  # tell lavaan that the repeat variable is binary
  fit <- sem(model, data = mydata, ordered="repeat.")
  
  # save the results for later statistical inference
  coefs <- parameterEstimates(fit)
  # paramsStand <- standardizedSolution(fit) # how are these standardized?
  
  # save
  data[s, 2] = coefs[['est']][1]
  data[s, 3] = coefs[['est']][2]
  data[s, 4] = coefs[['est']][3]
  data[s, 5] = coefs[['est']][12]
  data[s, 6] = coefs[['est']][13]
  data[s, 7] = coefs[['est']][14]
}

# put all together 
colnames(data) <- c('c', 'a', 'b', 'c1', 'indirect', 'direct', 'total')
write.csv(data, sprintf("%s/Data/CSV/SEMdata_lavaan.csv", mypath))


