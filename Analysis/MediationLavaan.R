# set path
mypath <- '/Users/anne/Data/pupilUncertainty_FigShare'

# activate the lavaan package
library("lavaan", lib.loc="~/Library/R/3.2/library")

# preallocate
nr_subjects = 27;
data = matrix(, nrow = nr_subjects, ncol = 9)

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
            decision_pupil ~ a1*uncertainty
            rt ~ a2*uncertainty

            # indirect effect through mediator
            repeat. ~ c1*uncertainty + cov2*stimrepeat + b1*decision_pupil + b2*rt

            # specify the effects we want to test
            indirectPup := a1*b1
            indirectRT  := a2*b2
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
  data[s, 6] = coefs[['est']][5]
  data[s, 7] = coefs[['est']][6]
  
  # also save the indirect effects directly
  data[s, 8] = coefs[['est']][16]
  data[s, 9] = coefs[['est']][17]

  # also test whether the indirect effect is significant for this person
  # data[s, 8] = coefs[['pvalue']][12]
  
  # also get some measure of the goodness of fit
  # gof = fitMeasures(fit)
  # data[s, 9] = gof['cfi']
  
}

# put all together 
colnames(data) <- c('c', 'a1', 'a2', 'c1', 'cov2', 'b1', 'b2', 'indirectPup', 'indirectRT')
write.csv(data, sprintf("%s/Data/CSV/SEMdata_lavaan.csv", mypath))


