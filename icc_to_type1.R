# load required packages
library(arm)

# set value
lvl2units <- 31
lvl1unitsperlvl2unit <- 40
targetICC <- .1 
sigma.a <- 1.7 ## use 2.5 for ICC 0.2; use 1.7 for 10
totalnumberoflines <- lvl2units * lvl1unitsperlvl2unit
subjectspercondition <- totalnumberoflines / 2 
ICCinfbound <- targetICC - 0.001
ICCsupbound <- targetICC + 0.001
j <- 0
counter <- 0
zr <- c(0, 0, 0, 0, 0, 0, 0)
rho <- 0.56
mu.a <- 0
mu.b <- 3
sigma.b <- 4
sigma.y <- 1
replications <- 100
group <- rep(1:lvl2units, rep(lvl1unitsperlvl2unit,lvl2units))
cond <- gl(2,subjectspercondition) 
lvl2groups <- gl(lvl2units,lvl1unitsperlvl2unit) 
Sigma.ab <- array(c(sigma.a^2, rho * sigma.a * sigma.b, rho * sigma.a * sigma.b, sigma.b^2), c(2, 2))
temp <- NULL
l <- 1
repeat { 
  ab <- mvrnorm(lvl1unitsperlvl2unit, c(mu.a, mu.b), Sigma.ab)
  a <- ab[,1]
  b <- ab[,2] 
  x <- rnorm(lvl2units*lvl1unitsperlvl2unit) 
  data <- rnorm(lvl2units*lvl1unitsperlvl2unit, a[group] + b*x,sigma.y)
  lm.fit <- lm(data~cond)
  sumary <- summary(lm.fit)
  cond.estimate <- sumary$coefficients[2] 
  cond.se <- sumary$coefficients[4]
  t.value <- cond.estimate / cond.se
  p.value <- sumary$coefficients[8] 
  lm.lvl2groups.fit <- lm(data~lvl2groups) 
  MSEwithin <- anova(lm.lvl2groups.fit)[2,2] / anova(lm.lvl2groups.fit)[2,1] 
  MSEgroups <- anova(lm.lvl2groups.fit)[1,2] / anova(lm.lvl2groups.fit)[1,1] 
  MSEbetween <- ((MSEgroups-MSEwithin) / lvl1unitsperlvl2unit) 
  ICC <- MSEbetween / (MSEwithin + MSEbetween)
  temp[l] <- ICC
  l <- l + 1
  t.Kish.correctn <- sqrt((1 + (lvl1unitsperlvl2unit-1) * ICC)) 
  Kishcorrctd.t.value <- t.value / t.Kish.correctn 
  Kishcorrectd.p.value <- 2 * pt(-abs(Kishcorrctd.t.value), anova(lm.fit)[2,1])
  if(ICC > ICCinfbound){ 
    if(ICC < ICCsupbound){ 
      counter <- counter + 1 
      zr <- c(zr, counter, t.value, p.value, ICC, t.Kish.correctn, Kishcorrctd.t.value, Kishcorrectd.p.value)}}
  j <- counter
  if (j > (replications - 1)) break } 

endofzr <- 7 * (replications + 1)
zrOK <- zr[8:endofzr]
tmp <- array(zrOK, c(7, counter))
t(tmp)

test <- t(tmp)

sum(test[, 3] < .05) / length(test[, 3]) * 100

mean(test[, 4])

    