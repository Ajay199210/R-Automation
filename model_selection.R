
# model_selection.R: A script that automates linear models algorithms selection.

#-- Load required libraries
library(boot)

#-------------------#
# Forward Selection #
#-------------------#

#-- Manual 'Forward Selection' on the nuclear data set with the 'cost' as the dep. var
head(nuclear)

nuc.0 <- lm(cost~1,data=nuclear)
summary(nuc.0)
add1(nuc.0,scope=.~.+date+t1+t2+cap+pr+ne+ct+bw+cum.n+pt,test="F")
nuc.1 <- update(nuc.0,formula=.~.+date)
summary(nuc.1)
add1(nuc.1,scope=.~.+date+t1+t2+cap+pr+ne+ct+bw+cum.n+pt,test="F")
nuc.2 <- update(nuc.1,formula=.~.+cap)
summary(nuc.2)
add1(nuc.2,scope=.~.+date+t1+t2+cap+pr+ne+ct+bw+cum.n+pt,test="F")
nuc.3 <- update(nuc.2,formula=.~.+pt)
summary(nuc.3)
add1(nuc.3,scope=.~.+date+t1+t2+cap+pr+ne+ct+bw+cum.n+pt,test="F")
nuc.4 <- update(nuc.3,formula=.~.+ne)
summary(nuc.4)
add1(nuc.4,scope=.~.+date+t1+t2+cap+pr+ne+ct+bw+cum.n+pt,test="F")
### cost ~ date + cap + pt + ne best fits the data in forward selection

#-- Manual 'Forward Selection' on the mtcars data set with the 'mpg' as the dep. var
head(mtcars)

car.0 <- lm(mpg~1,data=mtcar)
summary(car.0)
add1(car.0,scope=.~.+cyl+disp+hp+drat+wt+qsec+vs+am+gear+carb,test="F")
car.1 <- update(car.0,formula=.~.+wt)
summary(car.1)
add1(car.1,scope=.~.+cyl+disp+hp+drat+wt+qsec+vs+am+gear+carb,test="F")
car.2 <- update(car.1,formula=.~.+cyl)
summary(car.2)
add1(car.2,scope=.~.+cyl+disp+hp+drat+wt+qsec+vs+am+gear+carb,test="F")
### mpg ~ wt + cyl best fits the data in forward selection

#-- Automating the forward selection algorithm
fwd.select <- function(mydata, dep.var) {
  
  # Determining the dep & indep variables
  y <- which(colnames(mydata) == dep.var)
  x <- which(colnames(mydata) != dep.var) 
  
  # Re-arranging the columns with the dep var in the 1st col
  mydata <- mydata[, c(y,x)]
  
  # Building the formula scope for the forward selection using add1()
  vars = colnames(mydata)
  form.scope <- paste(paste(".~.", "+", paste(vars[2:length(vars)], collapse="+")))
  form.scope <- formula(form.scope)

  # Multiple models holder lists
  mul.mod.list <- list() # update()
  summary.list <- list() # summary(update)
  anova.df.list <- list() # add1()
  
  # Starting model formula for forward selection (overall intercept)
  mod0.form <- formula(mydata[,y] ~ 1)
  mul.mod.list[[1]] <- lm(mod0.form, data = mydata)
  
  # Applying the 'Forward Selection' algorithm
  for (i in 1:length(mydata)) {
    summary.list[[i]] <- summary(mul.mod.list[[i]])
    anova.df.list[[i]] <- add1(mul.mod.list[[i]], scope=form.scope, test="F")
    p.values <- anova.df.list[[i]]$`Pr(>F)`
    if (any(p.values < 0.05, na.rm = T)) {
      smallest.P <- min(p.values, na.rm=T)
      sig.var.ind <- which(anova.df.list[[i]]$`Pr(>F)` == smallest.P)
      # browser()
      sig.var.name <- rownames(anova.df.list[[i]])[sig.var.ind]
      new.form <- formula(paste(".~.", "+", sig.var.name, sep=" "))
      mul.mod.list[[i+1]] <- update(mul.mod.list[[i]], formula=new.form)
    } else {
      break
    }
  }
  return(summary.list)
}

#-- Testing the algorithm on the nuclear & mtcars data sets
nuc.fs <- fwd.select(nuclear,"cost")
nuc.fs

car.fs <- fwd.select(mtcars,"mpg")
car.fs

#--------------------------------#
# Backward Selection/Elimination #
#--------------------------------#

#-- Manual 'Backward Selection' on the nuclear data set with 'cost' as the dep. var
nuc.0 <- lm(cost~date+t1+t2+cap+pr+ne+ct+bw+cum.n+pt,data=nuclear)
summary(nuc.0)

drop1(nuc.0,test="F")
nuc.1 <- update(nuc.0,.~.-bw)
summary(nuc.1)

drop1(nuc.1,test="F")
nuc.2 <- update(nuc.1,.~.-pt)
summary(nuc.2)

drop1(nuc.2,test="F")
nuc.3 <- update(nuc.2,.~.-t1)
summary(nuc.3)

drop1(nuc.3,test="F")
nuc.4 <- update(nuc.3,.~.-ct)
summary(nuc.4)

drop1(nuc.4,test="F")
### cost ~ date + t2 + cap + pr + ne + cum.n best fits the data in backward elim

#-- Manual 'Backward Selection' on the mtcars data set with 'mpg' as the dep. var
max(coef(summary(car.7))[,4])

car.0 <- lm(mpg~cyl+disp+hp+drat+wt+qsec+vs+am+gear+carb,data=mtcars)
summary(car.1)

drop1(car.0,test="F")
car.1 <- update(car.0,.~.-cyl)
summary(car.1)

drop1(car.1,test="F")
car.2 <- update(car.1,.~.-vs)
summary(car.2)

drop1(car.2,test="F")
car.3 <- update(car.2,.~.-carb)
summary(car.3)

drop1(car.3,test="F")
car.4 <- update(car.3,.~.-gear)
summary(car.4)

drop1(car.4,test="F")
car.5 <- update(car.4,.~.-drat)
summary(car.5)

drop1(car.5,test="F")
car.6 <- update(car.5,.~.-drat)
summary(car.6)

drop1(car.6,test="F")
car.7 <- update(car.6,.~.-disp)
summary(car.7)

drop1(car.7,test="F")
car.8 <- update(car.7,.~.-hp)
summary(car.8)

drop1(car.8,test="F")
### mpg ~ wt + qsec + am best fits the data in backward selection

#-- Automating the backward selection algorithm
bwd.select <- function(mydata, dep.var) {
  
  # Determining the dep & indep variables
  y <- which(colnames(mydata) == dep.var)
  x <- which(colnames(mydata) != dep.var) 
  
  # Re-arranging the columns with the dep var in the 1st col
  mydata <- mydata[, c(y,x)]
  
  # Multiple models holder lists
  mul.mod.list <- list() # update()
  summary.list <- list() # summary(update)
  anova.df.list <- list() # drop1()
  
  # Starting model formula for backward elimination (fullest model)
  vars = colnames(mydata)
  mod0.form <- paste(paste(vars[1], "~"), paste(vars[2:length(vars)], collapse='+'))
  mod0.form <- formula(mod0.form)
  mul.mod.list[[1]] <- lm(mod0.form, data = mydata)
  
  # Applying the 'Backward Selection' algorithm
  for (i in 1:length(mydata)) {
    summary.list[[i]] <- summary(mul.mod.list[[i]])
    anova.df.list[[i]] <- drop1(mul.mod.list[[i]], test="F")
    p.values <- anova.df.list[[i]]$`Pr(>F)`
    if (any(p.values > 0.05, na.rm = T)) {
      largest.P <- max(p.values, na.rm=T)
      insig.var.ind <- which(anova.df.list[[i]]$`Pr(>F)` == largest.P)
      # browser()
      insig.var.name <- rownames(anova.df.list[[i]])[insig.var.ind]
      new.form <- formula(paste(".~.", "-", insig.var.name, sep=" "))
      mul.mod.list[[i+1]] <- update(mul.mod.list[[i]], formula=new.form)
    } else {
      break
    }
  }
  return(summary.list)
}

#-- Testing the algorithm on the nuclear & mtcars data sets
nuc.bw <- bwd.select(nuclear, "cost")
nuc.bw

car.bw <- bwd.select(mtcars, "mpg")
car.bw

#------------------------#
# Stepwise AIC Selection #
#------------------------#

#-- Manual 'AIC Selection' on the nuclear data set with 'cost' as the dep. var
car.null <- lm(mpg ~ 1, data = mtcars)
car.step <- step(car.null, scope=.~.+wt*hp*factor(cyl)*disp+am+
                   factor(gear)+drat+vs+qsec+carb)
summary(car.step)

#-- Automating the stepwise AIC selection algorithm
stepwise.aic <- function(mydata, dep.var, myscope) {
  
  # Determining the dependent variable (response/outcome)
  y <- which(colnames(mydata) == dep.var)
  
  # Building the null (starting) model
  depVar = colnames(mydata)[y]
  mod0.form <- formula(paste(depVar,"~",1))
  
  # Applying the 'Stepwise AIC Selection' algorithm 
  mydata.null <- lm(mod0.form, data = mydata)
  mydata.step <- step(mydata.null,myscope)
  
  return(summary(mydata.step))
}

#-- Testing the algorithm on the mtcars data set
car.aic <- stepwise.aic(mtcars, "mpg", .~.+wt*hp*factor(cyl)*disp+am+
                          factor(gear)+drat+vs+qsec+carb)

car.aic
### mpg ~ wt + hp + qsec + wt:hp bests fits the data according to the criterion
### based approach

