#---Final Project Presentation---
setwd("Desktop/UCSantaCruz/AMS 204/FinalProject")
sig_mm = read.table("mm6from11_signalProfile.txt", header = TRUE, sep="\t")
#str(sig_mm)
library(genefilter)
library(ggplot2)
library(viridis)


sig_mm$var<- rowVars(abs(sig_mm[c('pos1', 'pos2', 'pos3', 'pos4', 'pos5', 'pos6')]), na.rm=TRUE)
sig_mm$mean<- rowMeans(abs(sig_mm[,1:5]))
sig_mm$max <- apply(abs(sig_mm[c('pos1', 'pos2', 'pos3', 'pos4', 'pos5', 'pos6')]), 1, max)
sig.train = sig_mm[sample(nrow(sig_mm), 250000),]
sig.train[sig.train == -Inf] <- NA
sig.train <- sig.train[complete.cases(sig.train),]
sig.train1 = sig.train[sample(nrow(sig.train), 26666)]

sig.test = sig_mm[sample(nrow(sig_mm), 250000),]
sig.test$prop.log = log(sig.test$prop,10)
sig.test$prop.log[sig.test$prop.log == -Inf] <- -10000000000000000000000000000000000000000000
sig.test <- sig.test[complete.cases(sig.test),]
sig.s = sig_mm[(log(sig_mm[,7],10)>-4),]
sig.small = sig.s[(log(sig.s[,7],10)< -2),]
sig.small_train = sig.small_complete[sample(nrow(sig.small_complete), 250000),]
sig.small_test = sig.small[sample(nrow(sig.small_complete), 250000),]
sig.small_complete <- na.omit(sig.small)

sum(sig_mm$prop < 0.1)/nrow(sig_mm)
#summary(sig.train)
logit.prop = log(sig.train$prop/1-sig.train$prop)
logit = function(x){log(x/1-x)}
plot(sig.test$var,sig.test$prop,xlab="Variance",ylab="Proportion Miscalled 11-mers")
line(sig.test$fit.logM1, add = TRUE)
curve(predict(logM1,data.frame(var=x),type="response"),add=TRUE)
plot(logit(0:1), type='l',ylim = c(-4,4),xlim=c(0,1))
logit.prop.test = log(sig.test$prop)/log(1-sig.test$prop)
#----Get Density----
get_density <- function(x, y, n = 100) {
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

#-----Logistic Regression----
allM = glm(prop~., weight = total,family = binomial, data = sig.train)
step(allM,direction = "both")
step(allM, direction = "forward")
step(allM, direction = "backward")


logM = glm(sig.train$prop~sig.train$var, family = binomial, weight =sig.train$total )
summary(logM)
sig.test$fit = predict(logM,sig.test,type="response")
plot(log(sig.test$fit),log(sig.test$prop))
par(mfrow = c(2,1))
hist(log(sig.test$fit))
hist(log(sig.test$prop))
hist(sig.test$fit)
hist(sig.test$prop)

logM1 = glm(sig.train$prop~sig.train$var+sig.train$mean+sig.train$max, family = binomial, weight =sig.train$total )
summary(logM1)
step(logM1, direction = "both")
anova(logM1,test = "Chisq")
#AIC is a little bit better.

sig.test$fit.logM1 = predict(logM1,sig.test,type="response")
sig.test$logit.fit = predict(logM1,sig.test)
mean(sig.test$fit.logM1==sig.test$prop)
mean(sig.test$fit.logM1-sig.test$prop)/mean(sig.test$fit.logM1)



par(mfrow=c(3,1))
hist(sig.test$prop.log, breaks = 40, main = "", xlab = expression("log"[10]*"(Actual Proportion) with -Inf"))
hist(log(sig.test$prop,10), breaks = 40, main = "", xlab = expression("log"[10]*"(Actual Proportion) without -Inf"))
hist(log(sig.test$fit.logM1,10),main = "", xlim=c(-5,0),xlab = expression("log"[10]*"(Predicted Proportion)"))
mtext(expression(paste(bold("Frequencey of Actual and Model Proportions"))),cex=1.5, line = -3, outer = TRUE)

par(mfrow=c(1,1))
plot(density(log(sig.test$prop)),ylim=c(0,.4))
lines(density(log(sig.test$fit.logM1)),col=2)
plot(density(log(sig.test$fit.logM1)))

plot(log(sig.test$fit.logM1),(sig.test$prop))
#neither plot looks good
hist(log(sig.test$fit.logM1), main = "Estimated Logit", xlab = "log(estimated logit)")
hist(log(sig.test$logit), main = "Logit", xlab = "log(logit)")
hist(sig.test$logit.fit, main = "Estimated Logit", xlab = "Estimated Logit")
hist(sig.test$logit, main = "Logit", xlab = "Logit")
plot(log(sig.test$fit.logM1),log(sig.test$logit))
summary(sig.test$logit)
summary(sig.test$fit.logM1)
plot(sig.test$logit,sig.test$logit.fit)
plot(sig.test$prop, sig.test$logit.fit)
plot(sig.test$fit.logM1, sig.test$logit, main = "Est Logit vs Logit", xlab = "Est Logit", ylab = "Logit")
par(mfrow = c(1,3))
par(mfrow = c(1,1))
par(mfrow = c(1,2))
par(mfrow = c(2,1))

sig.test$logit = log(sig.test$prop/(1-sig.test$prop))


sig.test$density <- get_density(sig.test$logit.fit, sig.test$prop)

sig.test$densityM1 <- get_density(sig.test$fit.logM1, sig.test$prop)

ggplot(sig.test) + geom_point(aes(sig.test$fit.logM1, sig.test$prop))

ggplot(sig.test) + geom_point(aes(sig.test$logit,sig.test$prop))

ggplot(sig.test)+
  geom_point(aes(prop,logit.fit)) +
  geom_line(aes(prop,logit,col="Actual Logit"))+
  labs(title = "Predicted Logit Values for Mismatches, by Proportion Mismatched",
       x="Actual Proportion Mismatched",y="Logit")

ggplot(sig.test)+
  geom_point(aes(prop,fit.logM1))

#---log(prop)---
logM3 = glm(sig.train$log.prop~sig.train$var+sig.train$mean+sig.train$max, family = binomial, weight =sig.train$total )
sig.train$log.prop = log(sig.train$prop,10)
summary(sig.train$log.prop)
hist(sig.train$log.prop)
summary(logM1)
step(logM1, direction = "both")
anova(logM1,test = "Chisq")
#AIC is a little bit better..
sig.test$fit.logM1 = predict(logM1,sig.test,type="response")
sig.test$logit.fit = predict(logM1,sig.test)
mean(sig.test$fit.logM1==sig.test$prop)
mean(sig.test$fit.logM1-sig.test$prop)

#---data between 10^-4 and 10^-2---#
logM2 = glm(prop~var+mean+max, data = sig.small_train, family = binomial, weight =total )
summary(logM2)
step(logM2, direction = "both")
anova(logM2,test = "Chisq")
#AIC is a little bit better..
sig.small_test$fit.logM2 = predict(logM2,sig.small_test,type="response")

mean(sig.small_test$fit.logM2==sig.small_test$prop)
mean(sig.small_test$fit.logM2-sig.small_test$prop)/mean(sig.small_test$prop)

par(mfrow=c(2,1))
hist(log(sig.small_test$prop,10), breaks = 40, main = "", xlab = "Log Base 10 Transformed Actual Proportions")
hist(log(sig.small_test$fit.logM2,10),main = "", xlab = "Log Base 10 Transformed Predicted Proportions")
mtext(expression(paste(bold("Frequencey of Actual and Model Proportions"))),cex=1.5, line = -3, outer = TRUE)

par(mfrow=c(1,1))
plot(density(log(sig.small_test$prop)))
lines(density(log(sig.small_test$fit.logM2)),col=2)
plot(density(log(sig.small_test$fit.logM2)))

plot(log(sig.small_test$fit.logM2),(sig.small_test$prop))

#---Table---
library(xtable)
sig.test$fit.logM1 = predict(logM1,sig.test,type="response")

sig.test$cutprop = cut(sig.test$prop,
                     breaks=c(-1,0.001,0.006,0.02,1),
                     labels=c("~0","Low","Med.","High"))
sig.test$cutresponse = cut(sig.test$fit.logM1,
                         breaks=c(-1,0.001,0.006,0.02,1),
                         labels=c("M: ~0","Low","Med.","High"))
str(sig.test$cutresponse)
t = table(sig.test$cutprop,sig.test$cutresponse)
xtable(prop.table(t,margin=1),digits=c(3,3,3,3,3))

xtable(table(sig.test$cutprop,sig.test$cutresponse)) 


# #---change in signal over position---
delta.sig <- data.frame( header = c("pos1", "pos2",	"pos3",	"pos4",	"pos5",	"pos6"), wrong = c(67.880511,	69.364611,	84.882931,	102.731836,	86.035711,	74.961573), right = c(68.314212,	69.786145,	74.134501,	105.192991,	88.153018,	77.198632))
delta.sig$diff = abs(delta.sig$right-delta.sig$wrong)
str(delta.sig)
plot(delta.sig$header, delta.sig$diff,type = "l")
ggplot(data=delta.sig, aes(x=header, y=diff, group=1)) +
  geom_point()
str(sig.train)
alt.sig = data.frame(header = c("sig_diff","position"))
str(alt.sig)
#Actual Data#
library(tidyr)
sig.position <- sig_mm %>% gather(position, signal, pos1:pos6)
head(sig.position)
ggplot(sig.position, aes(x = position, y = abs(signal))) + geom_boxplot()

