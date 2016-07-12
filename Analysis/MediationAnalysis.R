# set path
mypath <- '/Users/anne/Data/pupilUncertainty_FigShare'

# activate the lavaan package
library("lavaan", lib.loc="~/Library/R/3.2/library")
library("mediation")

# example in the mediation package
set.seed(2014)
data("framing", package = "mediation")


med.fit <- lm(emo ~ treat + age + educ + gender + income, data = framing) # mediation model
out.fit <- glm(cong_mesg ~ emo + treat + age + educ + gender + income, data = framing, family = binomial("logit")) # outcome model

# fit the full  mediation model
med.out <- mediate(med.fit, out.fit, mediator = "emo", robustSE = TRUE, sims = 100)
summary(med.out)

for ( sj in 1 ) {
  # load data that have been coded to indicate repetition
  mydata = read.csv(sprintf("%s/Data/CSV/SEMdata_sj%02d.csv", mypath, sj))
  
  # make sure repeat is a logical array, so logistic regression is used
  mydata$repeat. <- factor(mydata$repeat.)
  
  # ============================================ #
  # MEDIATE PACKAGE
  # ============================================ #
  
  
  
  # ============================================ #
  # LAVAAN PACKAGE
  # ============================================ #

  # specify the mediation model
  model <- '# mediator
            decision_pupil ~ a*uncertainty
            repeat. ~ b*decision_pupil + c*uncertainty

            # specify the effects we want to test
            indirect := a*b
            direct   := c
            total    := c + (a*b)
            '
  
  # tell lavaan that the repeat variable is binary
  fit <- sem(model, data = mydata, ordered="repeat.")
  
  # display the results
  summary(fit)
  
  # save the results for later statistical inference
  params <- parameterEstimates(fit)
  paramsStand <- standardizedSolution(fit) # how are these standardized???
  
  # save
  write.csv(parans, sprintf("%s/Data/CSV/SEMdata_sj%02d.csv", mypath, sj))
}

