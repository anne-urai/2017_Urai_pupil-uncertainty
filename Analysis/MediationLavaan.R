# set path
mypath <- '/Users/anne/Data/pupilUncertainty_FigShare'

# activate the lavaan package
library("lavaan", lib.loc="~/Library/R/3.2/library")

# preallocate
nr_subjects = 27;
data = matrix(, nrow = nr_subjects, ncol = 6)

for ( s in 1:nr_subjects ) {
  
  # load data that have been coded to indicate repetition
  mydata = read.csv(sprintf("%s/Data/CSV/SEMdata_sj%02d.csv", mypath, s))
  
  # make sure repeat is a logical array, so logistic regression is used
  mydata$repeat. <- factor(mydata$repeat.)
  
  # specify the naive model
  model <- '# direct effect
			      repeat. ~ c*uncertainty + cov1*stimrepeat'
  
  fit <- sem(model, data=mydata, ordered=c("repeat."))
  coefs <- parameterEstimates(fit)
  gof = fitmeasures(fit)
  data[s, 1] = coefs[['est']][1]
  
  # specify the mediation model
  model <- '
            decision_pupil ~ a*uncertainty

            # indirect effect through mediator
            repeat. ~ c1*uncertainty + cov2*stimrepeat + b*decision_pupil

            # specify the effects we want to test
            indirect   := a*b
            # direct   := c1
            # total    := c1 + (a*b)
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
  data[s, 5] = coefs[['est']][4]
  
  # also save the indirect effects directly
  data[s, 6] = coefs[['est']][12]

}

# put all together 
colnames(data) <- c('c', 'a', 'c1', 'cov2', 'b','indirect')
write.csv(data, sprintf("%s/Data/CSV/SEMdata_lavaan.csv", mypath))


