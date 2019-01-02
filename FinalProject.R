setwd("Desktop/UCSantaCruz/AMS 204")
setwd("FinalProject")
library(ggplot2)
signal.data = read.table("template_median68pA.model.txt", header = TRUE)
ggplot(signal.data, aes(x = signal.data$level_mean)) + geom_density() + labs(x = "Current Level Mean")
summary(signal.data)

diff.data = read.table("hexamerSignalDiff.txt", sep = "\t", header = TRUE)
ggplot(diff.data, aes(x = diff.data$s_diff)) + geom_density() + labs(x = "Diff Mean Signal")

kmer = read.table("hexamer_sigdiff_nsig.txt", header = TRUE, sep="\t")
summary(kmer)

k_mm = data.frame(miscalled = kmer$miscalled.1.yes.0.no, sig_diff = abs(kmer$s_diff), n_sig = kmer$n_sig)
summary(k_mm)
k_mm.train = k_mm[sample(nrow(k_mm), 1000),]
k_mm.tst = k_mm[sample(nrow(k_mm), 1000),]
summary(k_mm.tst)
ggplot(k_mm, aes(x = sig_diff, y = miscalled)) + geom_boxplot()
(k_mm$sig_diff)
sapply(k_mm, function(x) sum(is.na(x)))
M = glm(formula = k_mm.train$miscalled ~ k_mm.train$n_sig, family = binomial(link='logit'))
M
summary(M)
fit = predict(M, k_mm.tst, type = "response")

pred_diff <- (fit - k_mm.tst$miscalled)

M2 = glm(formula = k_mm$miscalled ~ k_mm$n_sig + k_mm$sig_diff, family = binomial(link='logit'))
summary(M2)
which.max(k_mm[,2] )

M3 = glm(formula = k_mm$miscalled ~ 0 + k_mm$n_sig, family = binomial(link='logit'))
summary(M3)

M4 = glm(formula = k_mm$miscalled ~ k_mm$sig_diff, family = binomial(link='logit'))

#--------------------
el_mm = read.table("singleMismatchCount-11mer.txt", header = TRUE, sep="\t")
summary(el_mm)
el.df = data.frame(mis_kmer = el_mm$mis, ref_kmer = el_mm$ref11mer, mis_count = el_mm$misCount, cor_count = el_mm$corCount)
summary(el.df)
el.train = el.df[sample(nrow(el.df), 10000),]
el.tst = el.df[sample(nrow(el.df), 10000),]

hist(log(el.train$cor_count))
boxplot(el.train$mis_count)
plot(el.train$mis_count, el.train$cor_count, xlab = "Miscount", ylab = "Correct Count")

#----------------------------
base_mm = read.table("singleMismatchCount-11merByBase.txt", header = TRUE, sep="\t")

#----Linear Regression-------
sig_mm = read.table("mm6from11_signalProfile.txt", header = TRUE, sep="\t")
str(sig_mm)
library(genefilter)

sig_mm$var<- rowVars(sig_mm[c('pos1', 'pos2', 'pos3', 'pos4', 'pos5', 'pos6')], na.rm=TRUE)
sig_mm$mean<- rowMeans(abs(sig_mm[,1:5]))
sig_mm$max <- apply(sig_mm[c('pos1', 'pos2', 'pos3', 'pos4', 'pos5', 'pos6')], 1, max)
sig.train = sig_mm[sample(nrow(sig_mm), 250000),]
summary(sig.train)
logit.prop = log(sig.train$prop)/log(1-sig.train$prop)

M_L = lm(sig.train$prop~sig.train$var)
summary(M_L)
hist(M_L$residuals)
plot(M_L)

M_L1 = lm(sig.train$prop~sig.train$var+sig.train$mean)
plot(M_L1)
summary(M_L1)
hist(M_L1$residuals)

M_L2 = lm(sig.train$prop~sig.train$mean)
summary(M_L2)

M_L3 = lm(sig.train$prop~sig.train$max)
summary(M_L3)

M_L4 = lm(sig.train$prop~sig.train$var+sig.train$mean+sig.train$max)
summary(M_L4)
#max is not stat significant

M_L5 = lm(sig.train$prop~sig.train$mean+sig.train$max)
summary(M_L5)

M_L6 = lm(sig.train$prop~sig.train$var+sig.train$max)
summary(M_L6)


par(mfrow = c(3,2))
plot(M_L, which = 1)
plot(M_L, which = 2)
plot(M_L2, which = 1)
plot(M_L2, which = 2)
plot(M_L1, which = 1)
plot(M_L1, which = 2)

par(mfrow = c(3,2))
plot(M_L6, which = 1)
plot(M_L6, which = 2)
plot(M_L3, which = 1)
plot(M_L3, which = 2)
plot(M_L1, which = 1)
plot(M_L1, which = 2)

#----Linear Regression estimation of Logistic Regression (Data w/o 0's or 1's)--
zero_row_sub = apply(sig_mm, 1, function(row) all(row !=0 ))
sig_non0 = sig_mm[zero_row_sub,]
one_row_sub = apply(sig_non0, 1, function(row) all(row !=1 ))
sig_non01 = sig_non0[one_row_sub,]
#summary(sig_non01)
#str(sig_non01)
sig_non01$logit.prop = log(sig_non01$prop)/log(1-sig_non01$prop)

sig_non01.train = sig_non01[sample(nrow(sig_non01), 10000),]

M_non1 = lm(sig_non01.train$logit.prop~sig_non01.train$var)
summary(M_non1)
plot(M_non1)

M_non2 = lm(sig_non01.train$logit.prop~sig_non01.train$var+sig_non01.train$mean)
plot(M_non2)
summary(M_non2)

M_non3 = lm(sig_non01.train$logit.prop~sig_non01.train$mean)
summary(M_non3)
plot(M_non3)

#-----Logistic Regression----
sig.test = sig_mm[sample(nrow(sig_mm), 250000),]

logM = glm(sig.train$prop~sig.train$var, family = binomial, weight =sig.train$total )
summary(logM)
sig.test$fit = predict(logM,sig.test,type="response")
plot(log(sig.test$fit),log(sig.test$prop))
par(mfrow = c(2,1))
hist(log(sig.test$fit))
hist(log(sig.test$prop))
hist(sig.test$fit)
hist(sig.test$prop)

logM1 = glm(sig.train$prop~sig.train$var+sig.train$mean, family = binomial, weight =sig.train$total )
  lm(sig.train$prop~sig.train$var+sig.train$mean)
summary(logM1)
#AIC is a little bit better..
sig.test$fit.logM1 = predict(logM1,sig.test,type="response")
plot(sig.test$fit.logM1,sig.test$prop)
#neither plot looks good
step(logM,direction = "both")
step(logM, direction = "forward")
step(logM, direction = "backward")

#---Logistic Regressionw/o0---
sig_non0.train = sig_non0[sample(nrow(sig_non0), 10000),]
sig_non0.test = sig_non0[sample(nrow(sig_non0), 10000),]

logM.no0 = glm(sig_non0.train$prop~sig_non0.train$var, family = binomial, weight =sig_non0.train$total )
summary(logM)
sig_non0.test$fit = predict(logM,sig_non0.test,type="response")
plot(sig_non0.test$fit,sig_non0.test$prop)

logM1.no0 = glm(sig_non0.train$prop~sig_non0.train$var+sig_non0.train$mean, family = binomial, weight =sig_non0.train$total )
summary(logM1)
sig_non0.test$logM1.fit = predict(logM1,sig_non0.test,type="response")
plot(sig_non0.test$logM1.fit,sig_non0.test$prop)
