margsurv = function(indat, nreps.boot=100, bw.ini=0.5, bw.max=1, bw.srch=T, nreps.bw=30, grdwid=0.1, cnflvl=0.95, rnseed=43472) {

################################################################################################################################

# FUNCTION FOR MARGINAL SURVIVAL ESTIMATE FROM CASE-CONTROL FAMILY STUDY
# David Zucker and Malka Gorfine
#
# To use the function you have to load the following R packages:
#   Rcpp, survival
#
# Input arguments:
#
#   indat = data frame with one row for each relative and the following columns
#     fam.id = family id
#     prob.age = observation time of proband
#     prob.stat = proband status (1=case, 0=control)
#     rel.age = observation time of relative
#     rel.stat = relative status (1=event, 0=censored)
#
#   nreps.boot = number of outer bootstrap reps for variance estimation
#   bw.ini = initial bandwidth
#   bw.max = maximum bandwidth
#   bw.srch = flag for doing bandwidth search (TRUE or FALSE)
#   nreps.bw = number of inner bootstrap reps for bandwidth selection
#   grdwid = grid interval for bandwidth search
#   cnflvl = desired confidence level for percentile boostrap confidence interval
#
# Output: A list with the following three items:
#
#   out.table = a data frame with the following columns
#     TT = time
#     S_Hat = estimated marginal survival
#     SD = bootstrap estimate of standard deviation of S_Hat
#     SD_iqr = bootstrap estimate of interquartile range of S_Hat divided by the interquartile range of a N(0,1) r.v.
#     S_Lo = lower limit of percentile bootstrap confidence interval
#     S_Hi = upper limit of percentile bootstrap confidence interval
#
#   bw.fin = final selected bandwidth for original data
#
#   bw.fin.vec = vector of final selected bandwidths for the bootstrap samples
#
# Note: The marginal survival curve being estimated is the probability that a given individual in the population survives up
#   to time t given that his event time is greater than or equal to the lowest proband observation time

################################################################################################################################

### DATA SETUP ###

## SET UP CASE DATA
casedat0 = indat[which(indat$prob.stat==1),]
nrel = nrow(casedat0)
fam.id = casedat0$fam.id
idset = unique(fam.id)
nprob = length(idset)
prob.time = rep(0,nprob)
fam.id.new = rep(0,nrel)
nrel.vec = rep(0,nprob)
for (i in 1:nprob) {
  idx = which(fam.id==idset[i])
  prob.time[i] = casedat0$prob.age[idx[1]]
  fam.id.new[idx] = i
  nrel.vec[i] = length(idx)
}
prob.time.case = prob.time
nprob.case = nprob
nrel.case = nrel
nrel.vec.case = nrel.vec
reldat.case = cbind(fam.id.new,casedat0$prob.age,casedat0$rel.age,casedat0$rel.stat)

## SET UP CONTROL DATA
cntldat0 = indat[which(indat$prob.stat==0),]
nrel = nrow(cntldat0)
fam.id = cntldat0$fam.id
idset = unique(fam.id)
nprob = length(idset)
prob.time = rep(0,nprob)
fam.id.new = rep(0,nrel)
nrel.vec = rep(0,nprob)
for (i in 1:nprob) {
  idx = which(fam.id==idset[i])
  prob.time[i] = cntldat0$prob.age[idx[1]]
  fam.id.new[idx] = i
  nrel.vec[i] = length(idx)
}
prob.time.cntl = prob.time
nprob.cntl = nprob
nrel.cntl= nrel
nrel.vec.cntl = nrel.vec
reldat.cntl = cbind(fam.id.new,cntldat0$prob.age,cntldat0$rel.age,cntldat0$rel.stat)

TT = sort(unique(c(prob.time.case,prob.time.cntl)))
tmax = length(TT)

prob.time.case.new = rep(0,nprob.case)
prob.time.cntl.new = rep(0,nprob.cntl)
reldat.case.new = reldat.case
reldat.cntl.new = reldat.cntl

for (tt in 1:tmax) {
  tpcur = TT[tt]
  idx = which(prob.time.case==tpcur)
  prob.time.case.new[idx] = tt
  idx = which(reldat.case[,2]==tpcur)
  reldat.case.new[idx,2] = tt
  idx = which(prob.time.cntl==tpcur)
  prob.time.cntl.new[idx] = tt
  idx = which(reldat.cntl[,2]==tpcur)
  reldat.cntl.new[idx,2] = tt
}

prob.time.case = prob.time.case.new
prob.time.cntl = prob.time.cntl.new
reldat.case = reldat.case.new
reldat.cntl = reldat.cntl.new

### END OF DATA SETUP ###

#ADDTIONAL SETUPS RELATED TO BANDWIDTH SELECTION
bw.fin.vec = rep(0,nreps.boot)
bw.grid = seq(grdwid,bw.max,by=grdwid)
nbw.grid = length(bw.grid)
grhw = grdwid/2

#ADDITIONAL INITIAL SETUPS
S_Hat_Boot = matrix(0,nreps.boot,tmax)
iqrfac = -2*qnorm(0.25)
set.seed(rnseed)

## STORE ORIGINAL DATA ##
prob.time.cntl.orig = prob.time.cntl
reldat.cntl.orig = reldat.cntl
prob.time.case.orig = prob.time.case
reldat.case.orig = reldat.case

## RUN ESTIMATOR ON ORIGINAL DATA ##
estmte = est.proc(nprob.cntl, nrel.cntl, nrel.vec.cntl, prob.time.cntl, reldat.cntl,
  nprob.case, nrel.case, nrel.vec.case, prob.time.case, reldat.case, TT, tmax,
  bw.ini, bw.srch, bw.grid, nbw.grid, grhw, nreps.bw)
S_Hat = estmte$est.final
bw.fin = estmte$bw.fin

## BOOTSTRAP LOOP ##
for (ib in 1:nreps.boot) {
  prob.time.cntl.cur = NULL
  reldat.cntl.cur = NULL
  nrel.vec.cntl.cur = NULL
  prob.time.case.cur = NULL
  reldat.case.cur = NULL
  nrel.vec.case.cur = NULL
  index.cntl = sample.int(nprob.cntl, nprob.cntl, replace=T)
  for (iii in 1:nprob.cntl) {
    index.iii = index.cntl[iii]
    prob.time.cntl.cur = c(prob.time.cntl.cur, prob.time.cntl.orig[index.iii])
    rel.idx = which(reldat.cntl.orig[,1]==index.iii)
    nrel.cur = length(rel.idx)
    rel.iii = reldat.cntl.orig[rel.idx,2:4]
    if (nrel.cur==1) {
      rel.iii = c(iii,rel.iii)
    }
    else {
      rel.iii = cbind(rep(iii,nrel.cur),rel.iii)
    }
    reldat.cntl.cur = rbind(reldat.cntl.cur,rel.iii)
    nrel.vec.cntl.cur = c(nrel.vec.cntl.cur, nrel.cur)
  }
  index.case = sample.int(nprob.case, nprob.case, replace=T)
  for (iii in 1:nprob.case) {
    index.iii = index.case[iii]
    prob.time.case.cur = c(prob.time.case.cur, prob.time.case.orig[index.iii])
    rel.idx = which(reldat.case.orig[,1]==index.iii)
    nrel.cur = length(rel.idx)
    rel.iii = reldat.case.orig[rel.idx,2:4]
    if (nrel.cur==1) {
      rel.iii = c(iii,rel.iii)
    }
    else {
      rel.iii = cbind(rep(iii,nrel.cur),rel.iii)
    }
    reldat.case.cur = rbind(reldat.case.cur,rel.iii)
    nrel.vec.case.cur = c(nrel.vec.case.cur, nrel.cur)
  }
  nrel.cntl.cur = nrow(reldat.cntl.cur)
  nrel.case.cur = nrow(reldat.case.cur)
  estmte = est.proc(nprob.cntl, nrel.cntl.cur, nrel.vec.cntl.cur, prob.time.cntl.cur, reldat.cntl.cur,
    nprob.case, nrel.case.cur, nrel.vec.case.cur, prob.time.case.cur, reldat.case.cur, TT, tmax,
    bw.ini, bw.srch, bw.grid, nbw.grid, grhw, nreps.bw)
  S_Hat_Boot[ib,] = estmte$est.final
  bw.fin.vec[ib] = estmte$bw.fin
}

## WRAP-UP ##
ltail = (1-cnflvl)/2
rtail = cnflvl + ltail
SD = apply(S_Hat_Boot,2,sd,na.rm=T)
SD_iqr = apply(S_Hat_Boot,2,IQR,na.rm=T)/iqrfac
S_CI = apply(S_Hat_Boot,2,quantile,probs=c(ltail,rtail),type=8)
S_Lo = S_CI[1,]
S_Hi = S_CI[2,]
out.table = cbind(TT,S_Hat,SD,SD_iqr,S_Lo,S_Hi)
colnames(out.table) = c('TT','S_Hat','SD','SD_iqr','S_Lo','S_Hi')
out.table = as.data.frame(out.table)
ans = list(out.table=out.table, bw.fin=bw.fin, bw.fin.vec=bw.fin.vec)
return(ans)

}
### END OF MAIN FUNCTION margsurv ###

### AUXILIARY FUNCTIONS ########################################################################################################

# TRAPEZOID RULE INTEGRATION
my.integrate = function(xx,yy) {
  J0 = length(yy)-1
  xx1 = c(0,xx[1:J0])
  yy1 = c(0,yy[1:J0])
  base = xx-xx1
  height = (yy1+yy)/2
  ss = cumsum(base*height)
  #ss = c(0,ss)
  return(ss)
}

# FUNCTION TO COMPUTE Y AND dNR MATRICES
cppFunction('List Create_YdN(
  int nprob, int nrel, NumericMatrix reldat, NumericVector Unew, int lU) {

  // Yri and dNRi
  NumericMatrix Ymat(lU,nprob);
  NumericMatrix dNRmat(lU,nprob);
  for (int jj = 0; jj < lU; jj++) {
    double u = Unew(jj);
    for (int mm = 0; mm < nrel; mm++) {
      int ii = reldat(mm,0)-1;
      double obstim = reldat(mm,2);
      double status = reldat(mm,3);
        if (obstim >= u) {
          Ymat(jj,ii) += 1;
        }
        if ((obstim == u) && (status == 1)) {
          dNRmat(jj,ii) += 1;
        }
    }
  }

  return(Rcpp::List::create(
            Rcpp::Named("Ymat")=Ymat,
            Rcpp::Named("dNRmat")=dNRmat));

}')

# INTERCEPTS AND SLOPES FUNCTION
cppFunction('List Create_Alpha_Beta(
  int tmax, int nprob, NumericVector probtim, NumericMatrix Ymat, NumericMatrix dNRmat,
  NumericVector Te, NumericVector Unew, int lU, double width) {

  // Weights
  NumericMatrix Wmat(tmax,nprob);
  for (int ii = 0; ii < tmax; ii++) {
    double t = Te(ii);
    for (int rr = 0; rr < nprob; rr++) {
      Wmat(ii,rr) = 0;
      double z = (probtim(rr)- t)/width;
      if (fabs(z) <= 1) {
        Wmat(ii,rr) = pow((1-pow(z,2)),3);
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

    for (int rr = 0; rr < nprob; rr++) {
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
  int nprob, int lU, NumericMatrix Ymat, NumericMatrix dNRmat) {
  NumericVector surv(lU);
  double scur = 1;
  double s0, m0;
  for (int jj = 0; jj < lU; jj++) {
    s0 = 0;
    m0 = 0;
    for (int rr = 0; rr < nprob; rr++) {
      m0 += dNRmat(jj,rr);
      s0 += Ymat(jj,rr);
    }
    scur = scur * (1-m0/s0);
    surv(jj) = scur;
  }
  return(surv);
}')
## Combining Everything
Comp_Est = function(nprob.cntl, probtim.ctl, Ymat.ctl, dNRmat.ctl,
  nprob.case, probtim.case, Ymat.case, dNRmat.case, Te, Unew, width, tmax) {

  # Intercepts and Slopes
  lU = length(Unew)
  alpha.beta.case = Create_Alpha_Beta(tmax, nprob.case, probtim.case,
    Ymat.case, dNRmat.case, Te, Unew, lU, width)
  alpha.beta.control = Create_Alpha_Beta(tmax, nprob.cntl, probtim.ctl,
    Ymat.ctl, dNRmat.ctl, Te, Unew, lU, width)

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
bwsel = function(bw.grid, nprob.cntl, nrel.cntl, probtim.ctl, boot.reldat.cntl.arr,
   nprob.case, nrel.case, probtim.case, boot.reldat.case.arr, S.ini, Te, Unew, lU,
   tmax, nbw.grid, grhw, nreps.bw) {

  diff = array(rep(0,nreps.bw*nbw.grid*tmax),dim=c(nreps.bw,nbw.grid,tmax))
  umin = Unew[1]
  bw.max = bw.grid[nbw.grid]

  #LOOP OVER BOOTSTRAP DATASETS
  for (ir in 1:nreps.bw) {

    #GET BOOTSTRAP DATA
    boot.reldat.case = boot.reldat.case.arr[ir,1:nrel.case,1:4]
    boot.reldat.cntl = boot.reldat.cntl.arr[ir,1:nrel.cntl,1:4]

    #SET UP Y AND dNR MATRICES FOR CASES
    YdN = Create_YdN(nprob.case, nrel.case, boot.reldat.case, Unew, lU)
    Ymat.case = YdN$Ymat
    dNRmat.case = YdN$dNRmat

    #SET UP Y AND dNR MATRICES FOR CONTROLS
    YdN = Create_YdN(nprob.cntl, nrel.cntl, boot.reldat.cntl, Unew, lU)
    Ymat.ctl = YdN$Ymat
    dNRmat.ctl = YdN$dNRmat

    #KM Estimates
    km.cntl = KM(nprob.cntl, lU, Ymat.ctl, dNRmat.ctl)
    km.case = KM(nprob.case, lU, Ymat.case, dNRmat.case)
    km.cntl = c(rep(1,umin-1),km.cntl)
    km.case = c(rep(1,umin-1),km.case)
    km.cntl = km.cntl[1:tmax]
    km.case = km.case[1:tmax]

    #LOOP OVER BANDWIDTHS
    for (ibw in 1:nbw.grid) {
      bw.cur = bw.grid[ibw]
      est.cur.comp = Comp_Est(nprob.cntl, probtim.ctl, Ymat.ctl, dNRmat.ctl,
        nprob.case, probtim.case, Ymat.case, dNRmat.case, Te, Unew, bw.cur, tmax)
      est.cur = est.cur.comp$S
      est.cur = pmin(est.cur,km.cntl)
      est.cur = pmax(est.cur,km.case)
      diff.cur = est.cur - S.ini
      diff[ir,ibw,1:tmax] = diff.cur
    }

  }
  #END OF LOOP OVER BOOTSTRAP DATASETS

  #COMPUTE OBJECTIVE FUNCTION AND FIND BEST BW IN GRID
  bias = apply(diff,c(2,3),mean,na.rm=T)
  vari = apply(diff,c(2,3),var,na.rm=T)
  obj = bias^2 + vari
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
    boot.reldat.case = boot.reldat.case.arr[ir,1:nrel.case,1:4]
    boot.reldat.cntl = boot.reldat.cntl.arr[ir,1:nrel.cntl,1:4]

    #SET UP Y AND dNR MATRICES FOR CASES
    YdN = Create_YdN(nprob.case, nrel.case, boot.reldat.case, Unew, lU)
    Ymat.case = YdN$Ymat
    dNRmat.case = YdN$dNRmat

    #SET UP Y AND dNR MATRICES FOR CONTROLS
    YdN = Create_YdN(nprob.cntl, nrel.cntl, boot.reldat.cntl, Unew, lU)
    Ymat.ctl = YdN$Ymat
    dNRmat.ctl = YdN$dNRmat

    #KM Estimates
    km.cntl = KM(nprob.cntl, lU, Ymat.ctl, dNRmat.ctl)
    km.case = KM(nprob.case, lU, Ymat.case, dNRmat.case)
    km.cntl = c(rep(1,umin-1),km.cntl)
    km.case = c(rep(1,umin-1),km.case)
    km.cntl = km.cntl[1:tmax]
    km.case = km.case[1:tmax]

    #LOOP OVER THE NEW BW CANDIDATES
    for (ibw in 1:nbw.tmp) {
      bw.cur = bw.tmp[ibw]
      est.cur.comp = Comp_Est(nprob.cntl, probtim.ctl, Ymat.ctl, dNRmat.ctl,
        nprob.case, probtim.case, Ymat.case, dNRmat.case, Te, Unew, bw.cur, tmax)
      est.cur = est.cur.comp$S
      est.cur = pmin(est.cur,km.cntl)
      est.cur = pmax(est.cur,km.case)
      diff.cur = est.cur - S.ini
      diff.new[ir,ibw,1:tmax] = diff.cur
    }

  }

  #COMPUTE OBJECTIVE FUNCTION AND FIND BEST BW
  bias = apply(diff.new,c(2,3),mean,na.rm=T)
  vari = apply(diff.new,c(2,3),var,na.rm=T)
  obj = bias^2 + vari
  objbw.res.new = apply(obj,1,mean)
  bw.tmp = c(bw1,bw.tmp)
  objbw.tmp = c(objbw1,objbw.res.new)
  idxbw = which.min(objbw.tmp)
  bw.fin = bw.tmp[idxbw]
  if (idxbw == 1) {bias.fin = bias1}
  if (idxbw == 2) {bias.fin = bias[1,]}
  if (idxbw == 3) {bias.fin = bias[2,]}
  ans = list(bw.fin,bias.fin)
  names(ans) = c('bw.fin','bias.fin')

  return(ans)

}

## FULL ESTIMATION PROCEDURE ##
est.proc = function(nprob.cntl, nrel.cntl, nrel.vec.cntl, prob.time.cntl, reldat.cntl,
  nprob.case, nrel.case, nrel.vec.case, prob.time.case, reldat.case,
  TT, tmax, bw.ini, bw.srch, bw.grid, nbw.grid, grhw, nreps.bw) {

  ## LUMPED RELATIVES' SURVIVAL DATA ##
  rel.obs.times.case = reldat.case[,3]
  rel.status.case = reldat.case[,4]
  rel.obs.times.ctl = reldat.cntl[,3]
  rel.status.ctl = reldat.cntl[,4]
  rel.obs.times = c(rel.obs.times.case,rel.obs.times.ctl)
  rel.status = c(rel.status.case,rel.status.ctl)
  rel.ev.times = subset(rel.obs.times,rel.status==1)
  umin = min(rel.ev.times)
  umax = max(rel.ev.times)
  Unew = umin:umax
  lU = length(Unew)

  ## TIME TRANSFORMATION ##
  pbt.cntl = NULL
  for (i in 1:nprob.cntl) {
    pbt.cntl = c(pbt.cntl,rep(prob.time.cntl[i],nrel.vec.cntl[i]))
  }
  pbt.case = NULL
  for (i in 1:nprob.case) {
    pbt.case = c(pbt.case,rep(prob.time.case[i],nrel.vec.case[i]))
  }
  XPi1 = c(pbt.case,pbt.cntl)
  FF = ecdf(XPi1)
  ccc = sort(unique(XPi1))
  ccc = c(0,ccc)
  FFF = FF(ccc)
  Te = approx(ccc,FFF,1:tmax)$y
  prob.time.cntl.tr = Te[prob.time.cntl]
  prob.time.case.tr = Te[prob.time.case]

  ## Y AND dNR MATRICES ##
  YdN = Create_YdN(nprob.case, nrel.case, reldat.case, Unew, lU)
  probtim.case.orig = prob.time.case.tr
  Ymat.case.orig = YdN$Ymat
  dNRmat.case.orig = YdN$dNRmat
  YdN = Create_YdN(nprob.cntl, nrel.cntl, reldat.cntl, Unew, lU)
  probtim.ctl.orig = prob.time.cntl.tr
  Ymat.ctl.orig = YdN$Ymat
  dNRmat.ctl.orig = YdN$dNRmat

  #KM Estimates
  km.cntl = KM(nprob.cntl, lU, Ymat.ctl.orig, dNRmat.ctl.orig)
  km.case = KM(nprob.case, lU, Ymat.case.orig, dNRmat.case.orig)
  km.cntl = c(rep(1,umin-1),km.cntl)
  km.case = c(rep(1,umin-1),km.case)
  km.cntl = km.cntl[TT]/km.cntl[TT[1]-1]
  km.case = km.case[TT]/km.case[TT[1]-1]

  ## COMPUTATION OF INITIAL ESTIMATE ##
  est.ini = Comp_Est(nprob.cntl, probtim.ctl.orig, Ymat.ctl.orig, dNRmat.ctl.orig,
    nprob.case, probtim.case.orig, Ymat.case.orig, dNRmat.case.orig, Te, Unew, bw.ini, tmax)

  ## IF NO BANDWIDTH SEARCH, TAKE CURRENT ESTIMATE AS FINAL ##
  if (!bw.srch) {
    est.final = est.ini$S
    est.final = pmin(est.final,km.cntl)
    est.final = pmax(est.final,km.case)
    bw.fin = bw.ini
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

  bwboot = rep(0,nreps.bw)
  boot.reldat.cntl.arr = array(rep(0,nreps.bw*nrel.cntl*4),dim=c(nreps.bw,nrel.cntl,4))
  boot.reldat.case.arr = array(rep(0,nreps.bw*nrel.case*4),dim=c(nreps.bw,nrel.case,4))

  #GENERATE BOOTSTRAP SAMPLES
  for (ir in 1:nreps.bw) {

  boot.reldat.cntl = NULL
  boot.reldat.case = NULL

  for (i in 1:nprob.cntl) {
    ptm = prob.time.cntl[i]
    k = nrel.vec.cntl[i]
    boot.tim = sample.int(umax+1,k,replace=T,prob=p0[,ptm])
    boot.cens = sample.int(umax,k,replace=T,prob=pr.cens)
    boot.x = ifelse(boot.tim<=boot.cens, boot.tim, boot.cens)
    boot.d = ifelse(boot.tim<=boot.cens,1,0)
    for (ik in 1:k) {
      boot.reldat.cntl = rbind(boot.reldat.cntl,c(i,ptm,boot.x[ik],boot.d[ik]))
    }
  }

  for (i in 1:nprob.case) {
    ptm = prob.time.case[i]
    k = nrel.vec.case[i]
    boot.tim = sample.int(umax+1,k,replace=T,prob=p1[,ptm])
    boot.cens = sample.int(umax,k,replace=T,prob=pr.cens)
    boot.x = ifelse(boot.tim<=boot.cens, boot.tim, boot.cens)
    boot.d = ifelse(boot.tim<=boot.cens,1,0)
    for (ik in 1:k) {
      boot.reldat.case = rbind(boot.reldat.case,c(i,ptm,boot.x[ik],boot.d[ik]))
    }
  }

  boot.reldat.cntl.arr[ir,1:nrel.cntl,1:4] = boot.reldat.cntl
  boot.reldat.case.arr[ir,1:nrel.case,1:4] = boot.reldat.case

  }
  # END OF GENERATING BOOTSTRAP SAMPLES

  #BANDWIDTH SEARCH
  bw.sel.res = bwsel(bw.grid, nprob.cntl, nrel.cntl, probtim.ctl.orig, boot.reldat.cntl.arr,
    nprob.case, nrel.case, probtim.case.orig, boot.reldat.case.arr, S.ini, Te, Unew, lU,
    tmax, nbw.grid, grhw, nreps.bw)
  bw.fin = bw.sel.res$bw.fin

  # COMPUTE FINAL ESTIMATE
  est.fin = Comp_Est(nprob.cntl, probtim.ctl.orig, Ymat.ctl.orig, dNRmat.ctl.orig,
    nprob.case, probtim.case.orig, Ymat.case.orig, dNRmat.case.orig, Te, Unew, bw.fin, tmax)
  est.final = est.fin$S
  est.final = pmin(est.final,km.cntl)
  est.final = pmax(est.final,km.case)

  }
  ## END OF BANDWIDTH SEARCH PROCEDURE ##

  # RETURN RESULTS
  ans = list(est.final,bw.fin)
  names(ans) = c('est.final','bw.fin')
  return(ans)

}
## END OF FULL ESTIMATION PROCEDURE ##
