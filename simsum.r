#PROGRAM TO MERGE RESULTS FROM MARGINAL SURVIVAL SIMULATIONS
#AND PRODUCE SUMMARY STATS AND OUTPUT

library(stabledist)

runid = '31'
setwd("C:/Users/owner/Desktop/MS-results")

#SETTINGS - FIXED
n = 10^6
tmax = 110
iqrfac = -2*qnorm(0.25)
TT = 0:110

#SETTINGS - ADJUSTABLE
frailty = 'PS'      ## GAM=gamma, PS=positive stable
event.rate = 'LO'   ## end of study cuml event rate: HI=60%, LO=15%
tau = 1/3          ## within-family dependence parameter
wcut = 100          ## upper cutpoint for fraility (used only for PS)

#SURVIVAL PARAMETERS
mu = 0.01
pe = 4.6
cens.rate = 0.0055
if ((frailty=='GAM') & (event.rate=='HI')) {
  if (tau==1/3) {rate.par=1.0}
  if (tau==1/2) {rate.par=1.7}
}
if ((frailty=='GAM') & (event.rate=='LO')) {
  if (tau==1/3) {rate.par=0.115}
  if (tau==1/2) {rate.par=0.120}
}
if ((frailty=='PS') & (event.rate=='HI')) {
  if (tau==1/3) {rate.par=0.57}
  if (tau==1/2) {rate.par=0.53}
}
if ((frailty=='PS') & (event.rate=='LO')) {
  if (tau==1/3) {rate.par=0.041}
  if (tau==1/2) {rate.par=0.017}
}

if (frailty == 'GAM') {
  theta = 2*tau/(1-tau)
  gam = 1
}  
if (frailty == 'PS') {
  iota = complex(real=0, imaginary=1)
  theta = 1 - tau
  gam = abs(1-iota*tan(pi*theta/2))^(-1/theta)
} 

#TRUE MARGINAL SURVIVAL CURVE
if (frailty=='GAM') {
  Lam0 = rate.par * (mu*TT)^pe
  Lam.tru = log(1+theta*Lam0)/theta
  S_tq = exp(-Lam.tru) 
}
if (frailty=='PS') {
  w.ini1 = rstable(1.2*n,alpha=theta,beta=1,gamma=gam,delta=0,pm=1)
  w.ini2 = subset(w.ini1, w.ini1 <= wcut)
  w = w.ini2[1:n]
  time = (-log(1-runif(n,0,1))/(rate.par*w*(mu^pe)))^(1/pe)
  F.time = ecdf(time)
  S_tq = 1-F.time(TT)
  Lam.tru = -log(S_tq)
}  

#COMBINE DATA FROM DIFFERENT RUNS

wk.a = paste('Run',runid,'a.Rdata',sep='')
wk.b = paste('Run',runid,'b.Rdata',sep='')
wk.c = paste('Run',runid,'c.Rdata',sep='')
wk.d = paste('Run',runid,'d.Rdata',sep='')

e.a = new.env()
e.b = new.env()
e.c = new.env()
e.d = new.env()

load(wk.a, env=e.a)
load(wk.b, env=e.b)
load(wk.c, env=e.c)
load(wk.d, env=e.d)

S_Hat = rbind(
  e.a$S_Hat,
  e.b$S_Hat,
  e.c$S_Hat,
  e.d$S_Hat)
  
S_Hat_Boot_sd = rbind(
  e.a$S_Hat_Boot_sd,
  e.b$S_Hat_Boot_sd,
  e.c$S_Hat_Boot_sd,
  e.d$S_Hat_Boot_sd)
  
Lam_Hat = rbind(
  e.a$Lam_Hat,
  e.b$Lam_Hat,
  e.c$Lam_Hat,
  e.d$Lam_Hat)
  
Lam_Hat_Boot_sd = rbind(
  e.a$Lam_Hat_Boot_sd,
  e.b$Lam_Hat_Boot_sd,
  e.c$Lam_Hat_Boot_sd,
  e.d$Lam_Hat_Boot_sd) 
  
cover = rbind(
  e.a$cover,
  e.b$cover,
  e.c$cover,
  e.d$cover)         
  
bw.fin.vec = c(
  e.a$bw.fin.vec,
  e.b$bw.fin.vec,
  e.c$bw.fin.vec,
  e.d$bw.fin.vec) 

runs = length(bw.fin.vec)

## SUMMARY STATS ##

S_Hat_Mean = apply(S_Hat,2,mean)
S_Hat_Med = apply(S_Hat,2,median)
S_Mean_Bias = S_Hat_Mean - S_tq
S_Med_Bias = S_Hat_Med - S_tq
S_Hat_sd = apply(S_Hat,2,sd)
SD_Boot_Mean_S = apply(S_Hat_Boot_sd,2,mean)
SD_rat_S = SD_Boot_Mean_S/S_Hat_sd
S_Hat_sd_iqr = apply(S_Hat,2,IQR)/iqrfac

Lam_Hat_Mean = apply(Lam_Hat,2,mean)
Lam_Hat_Med = apply(Lam_Hat,2,median)
Lam_Mean_Bias = Lam_Hat_Mean - Lam.tru
Lam_Med_Bias = Lam_Hat_Med - Lam.tru
Lam_Hat_sd = apply(Lam_Hat,2,sd)
SD_Boot_Mean_Lam = apply(Lam_Hat_Boot_sd,2,mean)
SD_rat_Lam = SD_Boot_Mean_Lam/Lam_Hat_sd
Lam_Hat_sd_iqr = apply(Lam_Hat,2,IQR)/iqrfac

cover.rate = apply(cover,2,mean)

## SUMMARY OUTPUT ##

table = cbind(
  TT,
  S_tq,
  S_Hat_Mean,
  S_Hat_Med,
  S_Mean_Bias,
  S_Med_Bias,
  S_Hat_sd,
  SD_Boot_Mean_S,
  SD_rat_S,
  S_Hat_sd_iqr,
  cover.rate,
  Lam.tru,
  Lam_Hat_Mean,
  Lam_Hat_Med,
  Lam_Mean_Bias,
  Lam_Med_Bias,
  Lam_Hat_sd,
  SD_Boot_Mean_Lam,
  SD_rat_Lam,
  Lam_Hat_sd_iqr)

colnames(table) = c(
  "TT",
  "S_tq",
  "S_Hat_Mean",
  "S_Hat_Med",
  "S_Mean_Bias",
  "S_Med_Bias",
  "S_Hat_sd",
  "SD_Boot_Mean_S",
  "SD_rat_S",
  "S_Hat_sd_iqr",
  "Coverage",
  "Lam.tru",
  "Lam_Hat_Mean",
  "Lam_Hat_Med",
  "Lam_Mean_Bias",
  "Lam_Med_Bias",
  "Lam_Hat_sd",
  "SD_Boot_Mean_Lam",
  "SD_rat_Lam",
  "Lam_Hat_sd_iqr")

plot.file.name =  paste('Rplot_LL-',runid,'.pdf',sep='')
pdf(file=plot.file.name, onefile=T, paper='A4r')

plot(S_tq~TT, type="l", col=1, xlab= "Time", ylab = "Survival",
  main = paste("Marginal Survival Estimation", "\nrunid =", runid), ylim = c(0,1.2))
lines(S_Hat_Mean~TT, lty=2)
legend(x='bottomleft',y=NULL,legend = c("True", "Mean of Estimates"), lty = 1:2)

plot(S_tq~TT, type="l", col=1, xlab= "Time", ylab = "Survival",
  main = paste("Marginal Survival Estimation", "\nrunid =", runid), ylim = c(0,1.2))
lines(S_Hat_Med~TT, lty=2)
legend(x='bottomleft',y=NULL,legend = c("True", "Median of Estimates"), lty = 1:2)

hist(bw.fin.vec,freq=F)

dev.off()

tbl.file.name = paste('table_LL-',runid,'.csv',sep='')
write.csv(table,tbl.file.name)

print(summary(abs(S_Mean_Bias)))
print(summary(abs(S_Med_Bias)))
print(summary(S_Hat_sd))
print(summary(SD_Boot_Mean_S))
print(summary(SD_rat_S))
print(summary(S_Hat_sd_iqr))
print(summary(cover.rate))
print(summary(bw.fin.vec))
print(table(bw.fin.vec)/runs)

