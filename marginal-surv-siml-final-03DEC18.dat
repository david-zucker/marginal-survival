## PROGRAM FOR MARGINAL SURVIVAL ESTIMATE
## FROM CASE-CONTROL FAMILY STUDY
## SIMULATION CODE
## David Zucker and Malka Gorfine 
## December 2018

#setwd('C:/Users/owner/Desktop/MS-results')

#CLEAN UP WORKING ENVIRONMENT FROM PAST RUNS
rm(list = ls(), envir = .GlobalEnv)              #clear environment  
graphics.off()                                   #clear plots 
assign("last.warning", NULL, envir = baseenv())  #clear warnings
cat("\014")                                      #clear console

#LIBRARIES
library(Rcpp)
library(survival)
library(stabledist)
library(doParallel)
library(foreach)
library(doRNG)

## INITIAL SETUPS ##

wkfn = 'Run01a.Rdata'
batch = 'A'

# KEY SIMULATION PARAMETERS

frailty = 'GAM'      ## GAM=gamma, PS=positive stable
event.rate = 'HI'    ## end of study cuml event rate: HI=60%, LO=15%
n1 = 500             ## size of the sample
k = 1                ## number of relatives per proband
tau = 1/3            ## within-family dependence parameter
wcut = 100           ## upper cutpoint for fraility (used only for PS)

runs = 256           ## number of simulations   
nproc = 32           ## number of parallel machines

tmax = 110           ## maximum age
TT = 0:tmax          ## age range
n = 10^6             ## size of population from which the sample is taken

nreps.boot = 100     ## number of outer bootstrap reps for variance estimation

## BANDWIDTH SELECTION PARAMETERS ##
bw.ini = 0.5         ## initial bandwidth
bw.max = 1           ## maximum bandwidth 
bw.srch = T          ## flag for doing bandwidth search
nreps.bw = 30        ## number of inner bootstrap reps for bandwidth selection
grw = 0.1            ## grid interval for bandwidth search
nu = 1               ## tuning parameter for bias/variance tradeoff
bias_corr = F        ## flag for implementing the bootstrap bias correction 

# SURVIVAL DISTRIBUTION PARAMETERS

# mu and pe parameters as in Gorfine et al. (2016, Biostatistics)

# censoring is exponential with a rate of 0.55% per year 
# cumulative rate is 1 - exp(-110*0.0055) = 45%

# for HI scenario, rate.par is adjusted to yield an end-of-study survival of about 40%
# and a censoring rate (including end-of-study censoring) of about 60%

# for LO scenario, rate.par is adjusted to yield an end-of-study survival of about 85%

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

#ADDTIONAL SETUPS RELATED TO BANDWIDTH SELECTION
bw.grid = seq(grw,bw.max,by=grw)
nbw.grid = length(bw.grid)
grhw = grw/2

#SETUPS RELATED TO SIMULATIONS AND PARALLEL PROCESSING
if (batch=='A') {set.seed(912361)}
if (batch=='B') {set.seed(356998)}
if (batch=='C') {set.seed(309987)}
if (batch=='D') {set.seed(903984)}
doRNGversion("1.5.3")
if (nproc > 1) {registerDoParallel(nproc)}
runs.pp = runs/nproc

#ADDITIONAL INITIAL SETUPS
iqrfac = -2*qnorm(0.25)
tmax1 = tmax+1
S_Hat_Boot = matrix(0,nreps.boot,tmax1)
Lam_Hat_Boot = matrix(0,nreps.boot,tmax1)
icols = c((2:(k+1)),((k+3):(2*k+2)))
if (frailty == 'GAM') {
  theta = 2*tau/(1-tau)
  gam = 1
}  
if (frailty == 'PS') {
  iota = complex(real=0, imaginary=1)
  theta = 1 - tau
  gam = abs(1-iota*tan(pi*theta/2))^(-1/theta)
}  

###FUNCTIONS###

# TRAPEZOID RULE INTEGRATION
my.integrate = function(xx,yy) {
  J0 = length(yy)-1
  xx1 = c(0,xx[1:J0])
  yy1 = c(0,yy[1:J0])
  base = xx-xx1
  height = (yy1+yy)/2
  ss = cumsum(base*height)
  ss = c(0,ss)
  return(ss)
}

# FUNCTION TO COMPUTE Y AND dNR MATRICES
cppFunction('List Create_YdN(
  int n1, int k, NumericMatrix reldat, NumericVector Unew, int lU) {
  
  NumericMatrix Ymat(lU,n1);
  NumericMatrix dNRmat(lU,n1);
  for (int jj = 0; jj < lU; jj++) {
    double u = Unew(jj);
    for (int ii = 0; ii < n1; ii++) {
      Ymat(jj,ii) = 0;
      dNRmat(jj,ii) = 0;
      for (int ik = 0; ik < k; ik++) {
        if (reldat(ii,ik) >= u) {
          Ymat(jj,ii) += 1;
        }  
        if ((reldat(ii,ik) == u) && (reldat(ii,k+ik) == 1)) {
          dNRmat(jj,ii) += 1;
        }  
      }
    }
  }

  return(Rcpp::List::create(
            Rcpp::Named("Ymat")=Ymat,
            Rcpp::Named("dNRmat")=dNRmat));

}')

# INTERCEPTS AND SLOPES FUNCTION
cppFunction('List Create_Alpha_Beta(
  int tmax, int n1, int k, NumericVector probtim, NumericMatrix Ymat, NumericMatrix dNRmat,
  NumericVector kapwt, NumericVector Te, NumericVector Unew, int lU, double width) {
              
  // Weights
  NumericMatrix Wmat(tmax,n1);
  for (int ii = 0; ii < tmax; ii++) {
    double t = Te(ii);
    for (int rr = 0; rr < n1; rr++) {
      Wmat(ii,rr) = 0;
      double z = (probtim(rr)- t)/width;
      if (fabs(z) <= 1) {
        Wmat(ii,rr) = kapwt(rr)*pow((1-pow(z,2)),3);
      }     
    }
  } 
  
  // Computation of Simple Linear Regression

  NumericMatrix alf_mat_ini(lU,tmax);
  NumericMatrix beta_mat_ini(lU,tmax);
  double s0, s1, s2, m0, m1, t, xt, ss0, ss1, ss2, det;
  
  for (int ii = 0; ii < tmax; ii++) {
  for (int jj = 0; jj < lU; jj++) {
      
    t = Te(ii); 
    
    s0 = 0;
    s1 = 0;
    s2 = 0;
    m0 = 0;
    m1 = 0; 
    
    for (int rr = 0; rr < n1; rr++) {
      xt = probtim(rr)- t;
      m0 += Wmat(ii,rr)*dNRmat(jj,rr);
      m1 += Wmat(ii,rr)*dNRmat(jj,rr)*xt;
      ss0 = Wmat(ii,rr)*Ymat(jj,rr);
      ss1 = ss0*xt;
      ss2 = ss1*xt;
      s0 += ss0;
      s1 += ss1;
      s2 += ss2;
    }         
    
    alf_mat_ini(jj,ii) =  0;
    if (s0 > 0) {
      alf_mat_ini(jj,ii) =  m0/s0;
    }
    det = s0*s2 - pow(s1,2);
    if (det > 0) {
      alf_mat_ini(jj,ii) = (s2*m0-s1*m1)/det;
      beta_mat_ini(jj,ii) = (s0*m1-s1*m0)/det;
    }
    
  }
  }  
  
  // return results
  return(Rcpp::List::create(
    Rcpp::Named("alf.mat.ini")=alf_mat_ini,
    Rcpp::Named("beta.mat.ini")=beta_mat_ini));
  
}')

# KAPLAN-MEIER FUNCTION
cppFunction('NumericVector KM(
  int nprob, int lU1, NumericMatrix Ymat, NumericMatrix dNRmat,
  NumericVector kapwt) {
  NumericVector surv(lU1);
  double scur = 1;
  double s0, m0;
  for (int jj = 0; jj < lU1; jj++) {
    s0 = 0;
    m0 = 0;
    for (int rr = 0; rr < nprob; rr++) {
      m0 += kapwt(rr)*dNRmat(jj,rr);
      s0 += kapwt(rr)*Ymat(jj,rr);
    }
    scur = scur * (1-m0/s0);
    surv(jj) = scur;
  } 
  return(surv); 
}')

## Combining Everything 
Comp_Est = function(probtim.ctl, Ymat.ctl, dNRmat.ctl, kapwt.ctl, 
  probtim.case, Ymat.case, dNRmat.case, kapwt.case, Te, Unew, width) {

  # Intercepts and Slopes
  lU = length(Unew) 
  alpha.beta.case = Create_Alpha_Beta(tmax, n1, k, probtim.case, Ymat.case, dNRmat.case, kapwt.case, Te, Unew, lU, width) 
  alpha.beta.control = Create_Alpha_Beta(tmax, n1, k, probtim.ctl, Ymat.ctl, dNRmat.ctl, kapwt.ctl, Te, Unew, lU, width)
    
  # Integrand of Estimated Lambda 
  # rows are u and columns are t
  alf.mat.ini.case = alpha.beta.case$alf.mat.ini
  alf.mat.case = apply(alf.mat.ini.case, 2, cumsum)
  alf.mat.ini.ctl = alpha.beta.control$alf.mat.ini
  beta.mat.ini.ctl = alpha.beta.control$beta.mat.ini 
  alf.mat.ctl = apply(alf.mat.ini.ctl, 2, cumsum)
  beta.mat.ctl = apply(beta.mat.ini.ctl, 2, cumsum)
  S.mat.case.cum = exp(-alf.mat.case)  ## Conditional survival function for case
  S.mat.ctrl.cum = exp(-alf.mat.ctl)   ## Conditional survival function for control
  Lambda.cum.ctl.der = beta.mat.ctl    ## Derivative of control risk function 
  
  # Computation of Estimate
  numer = -S.mat.ctrl.cum * Lambda.cum.ctl.der
  denom = S.mat.ctrl.cum - S.mat.case.cum
  numer = numer*denom
  denom = denom^2
  numer.mean = apply(numer,2,mean)
  denom.mean = apply(denom,2,mean)
  lambda = rep(0,tmax)
  inc = which(denom.mean>0)
  lambda[inc] = numer.mean[inc]/denom.mean[inc]
  Lambda_Hat_tq = my.integrate(Te,lambda) 
  S_Hat = exp(-Lambda_Hat_tq)
  Lam = Lambda_Hat_tq
    
  # Return Results
  ans = list(S_Hat, Lam, S.mat.ctrl.cum, S.mat.case.cum)
  names(ans) = c('S','Lam','S0','S1')
  return(ans)

}

## BANDWIDTH SELECTION
bwsel = function(bw.grid, boot.ctl.sam.arr, boot.case.sam.arr, kapwt.ctl, kapwt.case, S.ini, Te, Unew, lU, lU1) {

  diff = array(rep(0,nreps.bw*nbw.grid*tmax),dim=c(nreps.bw,nbw.grid,tmax))
  umin = Unew[1]

  #LOOP OVER BOOTSTRAP DATASETS
  for (ir in 1:nreps.bw) {
    
    #GET BOOTSTRAP DATA
    boot.case.sample = boot.case.sam.arr[ir,1:n1,1:(2*k+2)]
    boot.control.sample = boot.ctl.sam.arr[ir,1:n1,1:(2*k+2)]
    
    #SET UP Y AND dNR MATRICES FOR CASES
    reldat = boot.case.sample[,icols]
    YdN = Create_YdN(n1, k, reldat, Unew, lU)
    probtim.case = boot.case.sample[,1]
    Ymat.case = YdN$Ymat
    dNRmat.case = YdN$dNRmat
    
    #SET UP Y AND dNR MATRICES FOR CONTROLS
    reldat = boot.control.sample[,icols]
    YdN = Create_YdN(n1, k, reldat, Unew, lU)
    probtim.ctl = boot.control.sample[,1]
    Ymat.ctl = YdN$Ymat
    dNRmat.ctl= YdN$dNRmat

    #KM Estimates
    km.cntl = KM(n1, lU1, Ymat.ctl, dNRmat.ctl, kapwt.ctl)
    km.case = KM(n1, lU1, Ymat.case, dNRmat.case, kapwt.case)
    km.cntl = c(rep(1,umin),km.cntl)
    km.case = c(rep(1,umin),km.case)
    
    #LOOP OVER BANDWIDTHS
    for (ibw in 1:nbw.grid) {
      bw.cur = bw.grid[ibw]
      est.cur.comp = Comp_Est(probtim.ctl, Ymat.ctl, dNRmat.ctl, kapwt.ctl, probtim.case, 
        Ymat.case, dNRmat.case, kapwt.case, Te, Unew, bw.cur)
      est.cur = est.cur.comp$S
      est.cur = pmin(est.cur,km.cntl)
      est.cur = pmax(est.cur,km.case)  
      diff.cur = est.cur - S.ini
      diff[ir,ibw,1:tmax] = diff.cur[2:tmax1]
    }
    
  }
  #END OF LOOP OVER BOOTSTRAP DATASETS

  #COMPUTE OBJECTIVE FUNCTION AND FIND BEST BW IN GRID  
  bias = apply(diff,c(2,3),mean,na.rm=T)
  vari = apply(diff,c(2,3),var,na.rm=T)
  obj = (bias^2)^nu + vari
  objbw.res = apply(obj,1,mean) 
  idxbw = which.min(objbw.res)   
  bw1 = bw.grid[idxbw]
  objbw1 = objbw.res[idxbw]
  bias1 = bias[idxbw,]
  
  #GET BW'S A HALF STEP DOWN AND A HALF STEP UP FROM CURRENT BEST
  bw2 = bw1 - grhw
  bw.tmp = bw2
  if (bw1 < bw.max) {
    bw3 = bw1 + grhw
    bw.tmp = c(bw2,bw3)
  }
  nbw.tmp = length(bw.tmp)
   
  diff.new = array(rep(0,nreps.bw*nbw.tmp*tmax),dim=c(nreps.bw,nbw.tmp,tmax))
   
  #LOOP OVER BOOTSTRAP DATASETS AGAIN
  for (ir in 1:nreps.bw) {
    
    #GET BOOTSTRAP DATA
    boot.case.sample = boot.case.sam.arr[ir,1:n1,1:(2*k+2)]
    boot.control.sample = boot.ctl.sam.arr[ir,1:n1,1:(2*k+2)]
    
    #SET UP Y AND dNR MATRICES FOR CASES
    reldat = boot.case.sample[,icols]
    YdN = Create_YdN(n1, k, reldat, Unew, lU)
    probtim.case = boot.case.sample[,1]
    Ymat.case = YdN$Ymat
    dNRmat.case = YdN$dNRmat
    
    #SET UP Y AND dNR MATRICES FOR CONTROLS
    reldat = boot.control.sample[,icols]
    YdN = Create_YdN(n1, k, reldat, Unew, lU)
    probtim.ctl = boot.control.sample[,1]
    Ymat.ctl = YdN$Ymat
    dNRmat.ctl= YdN$dNRmat

    #KM Estimates
    km.cntl = KM(n1, lU1, Ymat.ctl, dNRmat.ctl, kapwt.ctl)
    km.case = KM(n1, lU1, Ymat.case, dNRmat.case, kapwt.case)
    km.cntl = c(rep(1,umin),km.cntl)
    km.case = c(rep(1,umin),km.case)
    
    #LOOP OVER THE NEW BW CANDIDATES
    for (ibw in 1:nbw.tmp) {
      bw.cur = bw.tmp[ibw]
      est.cur.comp = Comp_Est(probtim.ctl, Ymat.ctl, dNRmat.ctl, kapwt.ctl, probtim.case,  
        Ymat.case, dNRmat.case, kapwt.case, Te, Unew, bw.cur)
      est.cur = est.cur.comp$S
      est.cur = pmin(est.cur,km.cntl)
      est.cur = pmax(est.cur,km.case)  
      diff.cur = est.cur - S.ini
      diff.new[ir,ibw,1:tmax] = diff.cur[2:tmax1]
    }
    
  }  
  
  #COMPUTE OBJECTIVE FUNCTION AND FIND BEST BW 
  bias = apply(diff.new,c(2,3),mean,na.rm=T)
  vari = apply(diff.new,c(2,3),var,na.rm=T)
  obj = (bias^2)^nu + vari
  objbw.res.new = apply(obj,1,mean) 
  bw.tmp = c(bw1,bw.tmp)
  objbw.tmp = c(objbw1,objbw.res.new)
  idxbw = which.min(objbw.tmp) 
  bw.fin = bw.tmp[idxbw]
  if (idxbw == 1) {bias.fin = bias1}
  if (idxbw == 2) {bias.fin = bias[1,]}
  if (idxbw == 3) {bias.fin = bias[2,]}
  bias.fin = c(0,bias.fin)
  ans = list(bw.fin,bias.fin)
  names(ans) = c('bw.fin','bias.fin') 
    
  return(ans) 

}

## FULL ESTIMATION PROCEDURE ##
est.proc = function(control.sample,case.sample,kapwt) {

  ## VECTORIZED VERSION OF RELATIVES' SURVIVAL DATA ## 
  rel.obs.times.case = as.vector(case.sample[,2:(k+1)])
  rel.status.case = as.vector(case.sample[,((k+3):(2*k+2))])
  rel.obs.times.ctl = as.vector(control.sample[,2:(k+1)])
  rel.status.ctl = as.vector(control.sample[,((k+3):(2*k+2))])
  rel.obs.times = c(rel.obs.times.case,rel.obs.times.ctl)
  rel.status = c(rel.status.case,rel.status.ctl)
  rel.ev.times = subset(rel.obs.times,rel.status==1)
  umin = min(rel.ev.times)
  umax = max(rel.ev.times)
  Unew = umin:umax
  lU = length(Unew)
 
  ## TIME TRANSFORMATION ##
  XPi = case.sample[,1]
  XPi1 = c(XPi,tmax)
  FF = ecdf(XPi1)
  ccc = sort(unique(XPi1))
  ccc = c(0,ccc)
  FFF = FF(ccc)
  Te = approx(ccc,FFF,1:tmax)$y
  case.sample.tr = case.sample
  control.sample.tr = control.sample
  case.sample.tr[,1] = Te[case.sample[,1]]
  control.sample.tr[,1] = Te[control.sample[,1]]
  
  ## Y AND dNR MATRICES ##
  reldat = case.sample.tr[,icols]
  YdN = Create_YdN(n1, k, reldat, Unew, lU)
  probtim.case.orig = case.sample.tr[,1]
  Ymat.case.orig = YdN$Ymat
  dNRmat.case.orig = YdN$dNRmat
  reldat = control.sample.tr[,icols]
  YdN = Create_YdN(n1, k, reldat, Unew, lU)
  probtim.ctl.orig = control.sample.tr[,1]
  Ymat.ctl.orig = YdN$Ymat
  dNRmat.ctl.orig = YdN$dNRmat

  #KM Estimates
  lU1 = length(umin:tmax)
  km.cntl = KM(n1, lU1, Ymat.ctl.orig, dNRmat.ctl.orig, kapwt)
  km.case = KM(n1, lU1, Ymat.case.orig, dNRmat.case.orig, kapwt)
  km.cntl = c(rep(1,umin),km.cntl)
  km.case = c(rep(1,umin),km.case)

  ## COMPUTATION OF INITIAL ESTIMATE ##
  est.ini = Comp_Est(probtim.ctl.orig, Ymat.ctl.orig, dNRmat.ctl.orig, kapwt, probtim.case.orig, 
    Ymat.case.orig, dNRmat.case.orig, kapwt, Te, Unew, bw.ini)
  
  ## IF NO BANDWIDTH SEARCH, TAKE CURRENT ESTIMATE AS FINAL ##
  if (!bw.srch) {
    est.final = est.ini$S
    bw.fin = bw.ini
    est.final = pmin(est.final,km.cntl)
    est.final = pmax(est.final,km.case)  
  }  
  
  ## EXECUTE BANDWIDTH SEARCH IF REQUESTED ##
  
  if (bw.srch) {

  # ESTIMATED S
  S.ini = est.ini$S
  S.ini = cummin(S.ini)
  S.ini = pmin(S.ini,km.cntl)
  S.ini = pmax(S.ini,km.case)  

  #PREPARE TO GENERATE BOOTSTRAP SAMPLES
  
  #ESTIMATE CENSORING DISTRIBUTION OF RELATIVES
  rel.status.cens = 1 - rel.status
  cens.surv = Surv(rel.obs.times,rel.status.cens)
  cens.km = survfit(cens.surv~1)
  kmsum = summary(cens.km)
  km.time = subset(kmsum$surv, kmsum$time <= umax)
  km.time = c(0,kmsum$time,umax+1)
  S.cens.ini = c(1,kmsum$surv,0)
  S.cens = approx(km.time,S.cens.ini,0:umax)$y
  S.cens[umax+1] = 0
  pr.cens = S.cens[1:umax] - S.cens[2:(umax+1)]
    
  # FURTHER PREPARATION FOR BOOTSTRAP SAMPLES

  # ESTIMATED S0 AND CORRESPONDING POINT PROBABILITIES   
  S0.ini = est.ini$S0
  S0.ini = rbind(rep(1,tmax),S0.ini)
  p0 = matrix(0,umax+1,tmax)
  for (it in 1:tmax) {
    S0.ini[,it] = cummin(S0.ini[,it])
    p0[Unew,it] = S0.ini[1:lU,it] - S0.ini[2:(lU+1),it]
    p0[umax+1,it] = S0.ini[lU+1,it]
  }  

  # ESTIMATED S1 AND CORRESPONDING POINT PROBABILITIES   
  S1.ini = est.ini$S1
  S1.ini = rbind(rep(1,tmax),S1.ini)
  p1 = matrix(0,umax+1,tmax)
  for (it in 1:tmax) {
    S1.ini[,it] = cummin(S1.ini[,it])
    p1[Unew,it] = S1.ini[1:lU,it] - S1.ini[2:(lU+1),it]
    p1[umax+1,it] = S1.ini[lU+1,it]
  }    

  # END PREPARE FOR BOOTSTRAP SAMPLES
  
  case.tim = case.sample[,1]
  bwboot = rep(0,nreps.bw)
  boot.ctl.sam.arr = array(rep(0,nreps.bw*n1*(2*k+2)),dim=c(nreps.bw,n1,(2*k+2)))
  boot.case.sam.arr = array(rep(0,nreps.bw*n1*(2*k+2)),dim=c(nreps.bw,n1,(2*k+2)))
   
  #GENERATE BOOTSTRAP SAMPLES
  for (ir in 1:nreps.bw) {
  
  boot.case.sample = matrix(0,n1,2*k+2)
  boot.case.sample[,1] = Te[case.tim]
  boot.case.sample[,k+2] = 1
  boot.control.sample = matrix(0,n1,2*k+2)
  boot.control.sample[,1] = Te[case.tim]
  boot.control.sample[,k+2] = 0

  for (i in 1:n1) {

    ptm = case.tim[i]
    
    boot.tim = sample.int(umax+1,k,replace=T,prob=p1[,ptm])
    boot.cens = sample.int(umax,k,replace=T,prob=pr.cens)
    boot.x = ifelse(boot.tim<=boot.cens, boot.tim, boot.cens)
    boot.d = ifelse(boot.tim<=boot.cens,1,0)
    boot.case.sample[i,2:(k+1)] = boot.x
    boot.case.sample[i,(k+3):(2*k+2)] = boot.d

    boot.tim = sample.int(umax+1,k,replace=T,prob=p0[,ptm])
    boot.cens = sample.int(umax,k,replace=T,prob=pr.cens)
    boot.x = ifelse(boot.tim<=boot.cens, boot.tim, boot.cens)
    boot.d = ifelse(boot.tim<=boot.cens,1,0)
    boot.control.sample[i,2:(k+1)] = boot.x
    boot.control.sample[i,(k+3):(2*k+2)] = boot.d
    
  } 
  
  boot.ctl.sam.arr[ir,1:n1,1:(2*k+2)] = boot.control.sample
  boot.case.sam.arr[ir,1:n1,1:(2*k+2)] = boot.case.sample
  
  }
  # END OF GENERATING BOOTSTRAP SAMPLES
  
  #BANDWIDTH SEARCH
  bw.sel.res = bwsel(bw.grid, boot.ctl.sam.arr, boot.case.sam.arr, kapwt, kapwt, S.ini, Te, Unew, lU, lU1)
  bw.fin = bw.sel.res$bw.fin
  
  # COMPUTE FINAL ESTIMATE
  est.fin = Comp_Est(probtim.ctl.orig, Ymat.ctl.orig, dNRmat.ctl.orig, kapwt, probtim.case.orig, 
      Ymat.case.orig, dNRmat.case.orig, kapwt, Te, Unew, bw.fin)
  est.final = est.fin$S
  est.final = pmin(est.final,km.cntl)
  est.final = pmax(est.final,km.case)  
  if (bias_corr) {est.final = est.final - bw.sel.res$bias.fin}

  }
  ## END OF BANDWIDTH SEARCH PROCEDURE ##

  # RETURN RESULTS
  ans = list(est.final,bw.fin)
  names(ans) = c('est.final','bw.fin')
  return(ans)

}
## END OF FULL ESTIMATION PROCEDURE ##

#DATA GENERATION FUNCTION
datgen = function(frailty,mu,pe,rate.par,cens.rate,theta,gam,wcut) {
  if (frailty == 'GAM') {w = rgamma(n,shape=(1/theta),scale=theta)}
  if (frailty == 'PS') {
    w.ini1 = rstable(1.2*n,alpha=theta,beta=1,gamma=gam,delta=0,pm=1)
    w.ini2 = subset(w.ini1, w.ini1 <= wcut)
    w = w.ini2[1:n]
  }   
  time = (-log(1-runif(n,0,1))/(rate.par*w*(mu^pe)))^(1/pe)
  for (ik in 1:k) {
    time = cbind(time,(-log(1-runif(n,0,1))/(rate.par*w*(mu^pe)))^(1/pe))
  }
  c = rexp(n,rate=cens.rate)
  for (ik in 1:k) {c=cbind(c,rexp(n,rate=cens.rate))}
  c = ifelse(c<tmax,c,tmax)  #censor at the maximum age
  x = ifelse(time<=c,time,c) #observed time
  d = ifelse(time<=c,1,0)    #event indicator
  x = ceiling(x) # for the sampling design of Davis
  y.dat = cbind(x,d) # x1 x2 x3 x4 d1 d2 d3 d4
  id.failpop  = seq(1:n)[y.dat[,k+2]==1]          #ids of case population
  id.case = sample(id.failpop, n1, replace=F)     #randomly select n1 cases
  case.sample = y.dat[id.case, ]                  #case proband families
  control.pop = subset(y.dat,y.dat[,k+2]==0)      #control population
  temp.age.con = control.pop[,1] #the following loop matches controls to cases by age of probands 
  id.match = vector()
  for (i in 1:n1) {
    id.match[i] = match(case.sample[i,1],temp.age.con)
    temp.age.con[id.match[i]] = NA
  }
  control.sample = control.pop[id.match,]

  F.time = ecdf(time)
  S_tq = 1-F.time(TT)

  ans = list(control.sample=control.sample, case.sample=case.sample, S_tq=S_tq)
  return(ans)

}
  
## SIMULATION FUNCTION ##
do.sim = function(runs.cur) {

S_Hat = matrix(0,runs.cur,tmax1)
Lam_Hat = matrix(0,runs.cur,tmax1)
S_Hat_Boot_sd  = matrix(0,runs.cur,tmax1)
Lam_Hat_Boot_sd  = matrix(0,runs.cur,tmax1)
cover = matrix(0,runs.cur,tmax1)
bw.fin.vec = rep(0,runs.cur)

ptmtm = proc.time()

for (r in 1:runs.cur) {
  
  ptn = proc.time()
  
  ## GENERATE DATA ##
  cur.dat = datgen(frailty,mu,pe,rate.par,cens.rate,theta,gam,wcut)
  control.sample = cur.dat$control.sample
  case.sample = cur.dat$case.sample

  ## STORE ORIGINAL DATA ##
  control.sample.orig = control.sample
  case.sample.orig = case.sample

  ## RUN ESTIMATOR ON ORIGINAL DATA ##
  kapwt = rep(1,n1)
  estmte = est.proc(control.sample,case.sample,kapwt)
  S_Hat_Cur = estmte$est.final 
 
  ## BOOTSTRAP LOOP ##
  for (ib in 1:nreps.boot) {
    idx.cur = sample.int(n1,n1,replace=T)
    control.sample.cur = control.sample.orig[idx.cur,]
    case.sample.cur = case.sample.orig[idx.cur,]
    estmte = est.proc(control.sample.cur,case.sample.cur,kapwt)
    S_Hat_Boot[ib,] = estmte$est.final
    Lam_Hat_Boot[ib,] = -log(estmte$est.final)
  }
  
  ## WRAP UP CURRENT REPLICATION ##
  Lam_Hat_Cur = -log(S_Hat_Cur)
  S_Hat[r,] = S_Hat_Cur
  Lam_Hat[r,] = Lam_Hat_Cur
  bw.fin.vec[r] = estmte$bw.fin
  SHBSD = apply(S_Hat_Boot,2,sd)
  LHBSD = apply(Lam_Hat_Boot,2,sd)
  S_Hat_Boot_sd[r,] = SHBSD
  Lam_Hat_Boot_sd[r,] = LHBSD
  S_Lo = apply(S_Hat_Boot,2,quantile,probs=0.025,type=8)
  S_Hi = apply(S_Hat_Boot,2,quantile,probs=0.975,type=8)
  cvr = (S_tq >= S_Lo) * (S_tq <= S_Hi)  
  cover[r,] = cvr 
  run.tim.sum = c(round(r),round(((proc.time() - ptmtm)[3])/60, digits = 2),
    round(((proc.time() - ptn)[3]), digits = 2))
  print(run.tim.sum)
  write(run.tim.sum, file='log.txt', append=T)
    
}
## END OF SIMULATION LOOP ##

ans = list(
  S_Hat=S_Hat,
  Lam_Hat=Lam_Hat,
  S_Hat_Boot_sd=S_Hat_Boot_sd,
  Lam_Hat_Boot_sd=Lam_Hat_Boot_sd,
  cover=cover,
  bw.fin.vec=bw.fin.vec)
return(ans)

}
## END OF SIMULATION FUNCTION ##

### MAIN PROGRAM ###

#TRUE MARGINAL SURVIVAL CURVE
if (frailty=='GAM') {
  Lam0 = rate.par * (mu*TT)^pe
  Lam.tru = log(1+theta*Lam0)/theta
  S_tq = exp(-Lam.tru) 
}
if (frailty=='PS') {
  simrun = datgen(frailty,mu,pe,rate.par,cens.rate,theta,gam,wcut)
  S_tq = simrun$S_tq
  Lam.tru = -log(S_tq)
}  

## RUN SIMULATIONS ##
my.sim = foreach(runs.cur=rep(runs.pp,nproc)) %dorng% do.sim(runs.cur)

## RETRIEVE AND COLLATE RESULTS ##
S_Hat = my.sim[[1]]$S_Hat
S_Hat_Boot_sd = my.sim[[1]]$S_Hat_Boot_sd
Lam_Hat = my.sim[[1]]$Lam_Hat
Lam_Hat_Boot_sd = my.sim[[1]]$Lam_Hat_Boot_sd
cover = my.sim[[1]]$cover
bw.fin.vec = my.sim[[1]]$bw.fin.vec
if (nproc>1) {
  for (cyc in 2:nproc) {
    S_Hat = rbind(S_Hat,my.sim[[cyc]]$S_Hat)
    S_Hat_Boot_sd = my.sim[[cyc]]$S_Hat_Boot_sd
    Lam_Hat = rbind(Lam_Hat,my.sim[[cyc]]$Lam_Hat)
    Lam_Hat_Boot_sd = rbind(Lam_Hat_Boot_sd,my.sim[[cyc]]$Lam_Hat_Boot_sd)
    cover = rbind(cover,my.sim[[cyc]]$cover)
    bw.fin.vec = c(bw.fin.vec,my.sim[[cyc]]$bw.fin.vec)
  }  
} 

#SAVE RESULTS
save.image(file=wkfn)
warnings() 
