library(corpcor); library(bayesSurv); library(MASS); library(GenKern); library(coda); library(psych)
library(mvtnorm); library(rootSolve); library(Matrix); library(MCMCpack); library(LaplacesDemon)
array_id = 1; set.seed(323);
array_id <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
iter = 100; burnin = 10; n = 100

closedform = TRUE
assumeXnormal = FALSE
difference = F
Xfit = T
condindfit = T
impute = T
method ='E'
SIM<-1500

alpha1 = 2; psi2 = 0
beta0 = 3; omega1 = 1
beta1 = 4.1; omega2 = 1
eps1 = 1; eps2 = 1; eps3 = 1; tau4 = 0.5
eta4 = 1

theta10<-0.15; theta11<-0.7; thetaT<-0.2142857
if(method == "B"){omega1 = 0; omega2 = 0}
if(method == "C"){omega1 = 1; omega2 = 0}
if(method == "D"){omega1 = 3; omega2 = 1; theta10 = 0.0766; theta11 = 0.3; thetaT = 0.255}
if(method == "E"){beta0 = 2; omega1 = -1; omega2 = 1; alpha1 = 3.8; psi2 = 1}
if(method == "F"){omega1 = -0.75; omega2 = 2}
#psi2 = 0
R<-matrix(rep(1,3*3),3,3); R[1,2] = R[2,1] = theta10;  R[1,3] = R[3,1] <-theta11;  R[2,3] = R[3,2] <- thetaT
R; eigen(R)
X = rnorm(n = iter*n, mean = eta4, sd = tau4)
X = rbinom(iter*n, 1, 0.5)

allsampmarg = mvrnorm(iter*n, mu = c(alpha1, beta0, beta1), R)
allsamp = allsampmarg + cbind(psi2 * X, omega1 * X, omega2 * X)

(theta11*eps3-theta10*eps2)/(eps1);
(beta1 - beta0)-((theta11*eps3-theta10*eps2)/(eps1))*(alpha1)
(beta1 + omega2 - beta0 - omega1)-((theta11*eps3-theta10*eps2)/(eps1))*(alpha1 + psi2)
(cor(allsamp[,1], allsamp[,3])*sd(allsamp[,3]) - cor(allsamp[,1], allsamp[,2])*sd(allsamp[,2]))/sd(allsamp[,1])
mean(allsamp[,3] - allsamp[,2]) - mean(allsamp[,1])*(cor(allsamp[,1], allsamp[,3])*sd(allsamp[,3]) - cor(allsamp[,1], allsamp[,2])*sd(allsamp[,2]))/sd(allsamp[,1])

P = matrix(rep(0,4*4),4,4); Q<-matrix(rep(1,3*3),3,3)

eta1 = alpha1
eta2 = beta0 + eta4 * omega1
eta3 = beta1 + eta4 * omega2
tau1 = eps1
tau2 = sqrt(eps2^2 + omega1^2 * tau4^2)
tau3 = sqrt(eps3^2 + omega2^2 * tau4^2)
pi10 = (eps2*theta10)/(sqrt(eps2^2+omega1^2*tau4^2))
pi11 = (eps3*theta11)/(sqrt(eps3^2+omega2^2*tau4^2))
piX0 = omega1*tau4/(eps2^2 + omega1^2 * tau4^2)
piX1 = omega2*tau4/(eps3^2 + omega2^2 * tau4^2)
#piX0 = cor(X, allsamp[,2])
#piX1 = cor(X, allsamp[,3])
#pi1X = cor(X, allsamp[,1])
pi1X = 0
piT = (omega1*omega2*tau4^2 * sqrt(eps2^2 + omega2^2 * tau4^2)* sqrt(
  eps3^2 + omega2^2 * tau4^2 ) + eps2 * eps3 * sqrt(eps2^2 + omega1^2 * tau4^2) *
    sqrt(eps3^2 + omega2^2 * tau4^2) * thetaT)/(eps2^2 * eps3^2 + eps3^2 * omega1^2 * 
   tau4^2 + eps2^2 * omega2^2 * tau4^2 + omega1^2 * omega2^2 * tau4^4)

Q[1,2] = Q[2,1] = pi10;  Q[3,1] = Q[1,3]<-pi11; Q[3,2] = Q[2,3]<-piT

thetaT - theta10/theta11

P[1,1] = tau1^2; P[2,2] = tau2^2; P[3,3] = tau3^2; P[4,4] = tau4^2
P[1,2] = P[2,1] = tau1*tau2*pi10; P[1,3] = P[3,1] = tau1*tau3*pi11
P[1,4] = P[4,1] = tau1*tau4*pi1X; P[2,3] = P[3,2] = tau2*tau3*piT
P[2,4] = P[4,2] = tau2*tau4*piX0; P[3,4] = P[4,3] = tau3*tau4*piX1
P; eigen(P)

round(ginv(P[1:3, 1:3]), 5); piT - pi10/pi11
V = diag(c(tau1, tau2, tau3))
(pi11*tau3-pi10*tau2)/(tau1); (eta3 - eta2)-((pi11*tau3-pi10*tau2)/(tau1))*(eta1)

S<-diag(c(eps1, eps2, eps3)); 

allsampmarg = mvrnorm(iter*n, mu = c(eta1, eta2, eta3, eta4), P)
trueR = R

P = diag(sqrt(diag(ginv(R))))
Sigdiag = P^2; Sig<- S%*%R%*%S; C = ginv(P%*%Sig%*%P)
Siginv<-P%*%C%*%P; SIG<-S%*%R%*%S; SIGinv<-P%*%C%*%P
a = b = 0.1

lm(allsamp[,1] ~ X); lm(allsamp[,2] ~ X); lm(allsamp[, 3] ~ X)

holdmu<-matrix(rep(0,3*SIM),3,SIM); holdmu1<-matrix(rep(0,3*SIM),3,SIM)
holdS<-array(rep(0,3*3*SIM),dim=c(3,3,SIM)); holdT<-array(rep(0,3*3*SIM),dim=c(3,3,SIM))
holdR<-array(rep(0,3*3*SIM),dim=c(3,3,SIM)); holdR[,,1]=R; holdS[,,1]=S; holdT[,,1]=ginv(S)
holdC<-array(rep(0,3*3*SIM),dim=c(3,3,SIM)); holdpsi1<-(rep(0,1*SIM)); holdpsi2<-(rep(0,1*SIM))
holdomega1<-(rep(0,1*SIM)); holdomega2<-(rep(0,1*SIM))
holdalpha0<-(rep(0,1*SIM)); holdalpha1<-(rep(0,1*SIM)); 
holdbeta0<-(rep(0,1*SIM)); holdbeta1<-(rep(0,1*SIM)); 
Sigmacoef=matrix((NA),SIM*6,nrow=6); holdomega1[1] = omega1; holdomega2[1] = omega2; holdpsi2[1] = psi2

if(difference == TRUE){  
  allsampd = allsampmarg
   allsampd = allsamp
 # allsampd[,1:3] = allsampd[,1:3] + cbind(psi2 * X, omega1 * X, omega2 * X)
 allsampd[,2:3] = allsampd[, 2:3] -  X
  
  sigs1 = sd(allsampd[,1])
  sigt0 = sqrt(tau2^2 + tau4^2 -2*tau2*tau4*piX0)
  sigt1 = sqrt(tau3^2 + tau4^2 -2*tau3*tau4*piX1)

  #sigt0 = sd(allsampd[,2])
  #sigt1 = sd(allsampd[,3])
  rho10 = cor(allsampd[,1], allsampd[,2])
  rho11 = cor(allsampd[,1], allsampd[,3])
  rhot = cor(allsampd[,2], allsampd[,3])
  (tau1*tau2*pi10 - tau1*tau4*pi1X)/(sigs1*sigt0)
  (tau1*tau3*pi11 - tau1*tau4*pi1X)/(sigs1*sigt0)
  
  (tau2*tau3*piT - tau2*tau4*piX0 - tau3*tau4*piX1 + tau4^2)/(sigt0*sigt1)
  cor(allsampd[,2], allsampd[,3])
  
  alpha1m = mean(allsampd[,1])
  beta0m = mean(allsampd[,2])
  beta1m = mean(allsampd[,3])
  eta2 - eta4
  eta3 - eta4
  rhot
  rho10/rho11
  Q[1,2] = Q[2,1] = rho10;  Q[1,3] = Q[3,1] = rho11;  Q[2,3] = Q[3,2] = rhot

  (theta11*eps3-theta10*eps2)/(eps1);

  (rho11*sigt1 - rho10*sigt0)/sigs1
  (beta1m - beta0m) - alpha1m*((rho11*sigt1 - rho10*sigt0)/sigs1)

  S<-diag(c(sigs1, sigt0, sigt1)); 
  holdS[,,1] = diag(c(sigs1, sigt0, sigt1))
  holdR[,,1] = Q
  holdmu[,1] = c(alpha1m, beta0m, beta1m)
  allsamp = allsampd
  }

alltrt<-c(rep(0,iter*n)); for (t in 1:iter){  alltrt[(t*n-(n/2-1)):(n*t)]<-1}
setind <- as.vector(mapply(rep, 1:iter, n)); set <- as.numeric(array_id);  X = X[setind == set]
samp=allsamp[setind==set,]; trt=alltrt[setind==set]
lm(samp[,1] ~ X); lm(samp[,2] ~ X); lm(samp[,3] ~ X)

sobs<-c(samp[(n/2):n,2]); tobs<-c(samp[1:(n/2-1),1],samp[(n/2):n,2])

ST<-samp; ST[1:(n/2),1]<-0; ST[(n/2+1):n,2]<-0; ST[1:(n/2),3]<-0; save=ST

fsig<-function(n,s,j){
  return(-(n)*log(s)+(-0.5*j)			)}

fdelt<-function(n, R, j){
  return(-(n/2)*log(det(R))+(-0.5*j)				)}

fdeltBeta<-function(n, R, j){
  p = 5; q = 6 # want p>q
  a = -0.4 # lower bound of beta
  b = 1 # upper bound of beta
  return( -(n/2) * log(det(R)) + (-0.5 * j) + (q - 1) * log(b - R[2,3]) + # beta
            (p - 1) * log(R[2,3] - a) ) # alpha
}

meanBeta = (5*1 + 6*(-0.4))/(5+6)

s1=diag(S)[1]; s2=diag(S)[2]; s3=diag(S)[3]
SIG0<-diag(rep(10^6,3)); SIG0S<-rep(.01,2)
SIG0Smat<-diag(.01,1)


if(Xfit == TRUE){
  holdmu[,1] = c(alpha1, beta0, beta1)
  holdmu1[,1] = c(alpha1 + psi2, beta0 + omega1, beta1 + omega2)
}


if(difference == FALSE & Xfit == FALSE){  holdmu[,1] = c(eta1, eta2, eta3)}

S1 = seq(min(ST[, 1]), max(ST[,1]), 0.1)

prob_t1t0<-matrix(0,(SIM), length(S1))

tauSS0=tauSS1=tauST0=tauST1=c(1); SIGSS0=SIGSS1=SIGST0=SIGST1=diag(1);
SIG0S<-diag(rep(10^6,2))
v0=3; V0=1
Xmat=cbind(1,X)
sim=2

surrogprob = function(x, S1, gamma0, gamma1, rho10, rhot, rho11, sigs1, sigt0, sigt1){
  pnorm(x, gamma0 + gamma1*S1,  sd = sqrt(sigt0^2 - sigt0^2*rho10^2 + sigt1^2 - 
                                            sigt1^2*rho11^2 - 2*rhot*sigt0*sigt1 + 2*rho10*rho11*sigt0*sigt1), lower.tail = FALSE)}

rBeta_ab <- function(n, shape1=2, shape2=3, a = 0, b = 1, params = list(shape1, shape2, a, b),...){
  if(!missing(params)){
    shape1 <- params$shape1;  shape2 <- params$shape2
    a <- params$a;   b <- params$b
  }
  X <- rbeta(n,shape1,shape2);  out <- (b-a)*X + a
  return(out)
}

slope=intmarginalizeclosed=slopemarginalizeclosed=slopemarginalize=array(0,c((sim-1),1))
slopemarginalizewtcl = intmarginalizewtcl = intmarginalize= (rep(0,1*SIM)); 
int= int1 =array(0,c((sim-1),1))

slope[1] = (holdR[1,3,1]*holdS[3,3,1]-holdR[1,2,1]*holdS[2,2,1])/(holdS[1,1,1])
(holdmu[3,1]-holdmu[2,1])-slope[1]*(holdmu[1,1])
(holdmu1[3,1]-holdmu1[2,1])-slope[1]*(holdmu1[1,1])
slope[1]

treateffects <- data.frame(margFX = numeric(1), margFXSE = numeric(1),
                           condFX = numeric(1), condFXSE = numeric(1)      
)

treateffects[1] = coef(lm(c(ST[trt==0,2], ST[trt==1,3]) ~ trt))[2]
treateffects[2] = summary(lm(c(ST[trt==0,2], ST[trt==1,3]) ~ trt))$coef[2,2]
treateffects[3] = coef(lm(c(ST[trt==0,2], ST[trt==1,3]) ~ trt + c(X[trt==0], X[trt==1])))[2]
treateffects[4] = summary(lm(c(ST[trt==0,2], ST[trt==1,3]) ~ trt + c(X[trt==0], X[trt==1])))$coef[2,2]

treateffects[1] = coef(lm((ST[c(TRUE, FALSE), 2] -  ST[c(TRUE, FALSE), 3]) ~ ST[c(TRUE, FALSE), 1]))[2]
treateffects[2] = summary(lm((ST[c(TRUE, FALSE), 2] -  ST[c(TRUE, FALSE), 3]) ~ ST[c(TRUE, FALSE), 1]))$coef[2,2]
treateffects[3] = coef(lm((ST[c(TRUE, FALSE), 2] -  ST[c(TRUE, FALSE), 3]) ~ ST[c(TRUE, FALSE), 1] + X[c(T, F)]))[2]
treateffects[4] = summary(lm((ST[c(TRUE, FALSE), 2] -  ST[c(TRUE, FALSE), 3]) ~ ST[c(TRUE, FALSE), 1] + X[c(T, F)]))$coef[2,2]

treateffects <- data.frame(margFX = numeric(1), margFXSE = numeric(1),
                           slope = numeric(1), condFXSE = numeric(1) ,
                           int1 = numeric(1), int2 = numeric(1),
                           slope1 = numeric(1), slope2 = numeric(1)
)

treateffects[1] = coef(lm(c(ST[trt==0,2], ST[trt==1,3]) ~ trt))[2]
treateffects[2] = summary(lm(c(ST[trt==0,2], ST[trt==1,3]) ~ trt))$coef[2,2]
treateffects[3] = coef(lm(c(ST[trt==0,2], ST[trt==1,3]) ~ trt + c(X[trt==0], X[trt==1])))[2]
treateffects[4] = summary(lm(c(ST[trt==0,2], ST[trt==1,3]) ~ trt + c(X[trt==0], X[trt==1])))$coef[2,2]

treateffects[5] = coef(lm((ST[c(TRUE, FALSE), 2] -  ST[c(TRUE, FALSE), 3]) ~ ST[c(TRUE, FALSE), 1]))[1]
treateffects[6] = summary(lm((ST[c(TRUE, FALSE), 2] -  ST[c(TRUE, FALSE), 3]) ~ ST[c(TRUE, FALSE), 1]))$coef[1,2]
treateffects[7] = summary(lm((ST[c(TRUE, FALSE), 2] -  ST[c(TRUE, FALSE), 3]) ~ ST[c(TRUE, FALSE), 1]))$coef[2]
treateffects[8] = summary(lm((ST[c(TRUE, FALSE), 2] -  ST[c(TRUE, FALSE), 3]) ~ ST[c(TRUE, FALSE), 1]))$coef[2,2]

fname <- paste('covs','treateffects', n, method, array_id,'.txt',sep="")
write.table(treateffects, file=fname, sep="\t",  col.names=T)
#break

lm(samp[,1] ~ X); lm(samp[,2] ~ X); lm(samp[, 3] ~ X)

while(sim<=SIM){
  mu=cbind(holdalpha1[sim-1] + holdpsi2[sim-1]*X, holdbeta0[sim-1] + holdomega1[sim-1]*X,  holdbeta1[sim-1] + holdomega2[sim-1]*X)
  #mu=cbind(holdalpha1[sim-1], holdbeta0[sim-1],  holdbeta1[sim-1])
  Sig = holdS[,,sim-1]%*%holdR[,,sim-1]%*%holdS[,,sim-1]
if(impute == TRUE){
      
  for(i in 1:n){
    if(trt[i]==0){
      ST[i,c(1,3)]<-c(mu[i,c(1,3)]+(Sig[c(1,3),c(2)]%*%ginv(Sig[c(2),c(2)]))%*%(ST[i,c(2)]-mu[i,c(2)]))+
        mvrnorm(1,c(0,0),Sig[c(1,3),c(1,3)]-Sig[c(1,3),c(2)]%*%ginv(Sig[c(2),c(2)])%*%Sig[c(2),c(1,3)])
    }
    if(trt[i]==1){
      ST[i,c(2)]<-c(mu[i,c(2)]+(Sig[c(2),c(1,3)]%*%ginv(Sig[c(1,3),c(1,3)]))%*%(ST[i,c(1,3)]-mu[i,c(1,3)]))+
        mvrnorm(1, c(0),Sig[c(2),c(2)]-Sig[c(2),c(1,3)]%*%ginv(Sig[c(1,3),c(1,3)])%*%Sig[c(1,3),c(2)])
    }
  }
  
  if(sim == 2){ST = samp}
  
  #### Estimate S1, T0, T1|X
  Xmat = rep(1, n) # n by k
  Sigmats = Sig # 3 x 3

  Lambda0t = matrix(c(rep(0, 1)), nrow = 1); diag(Lambda0t) = .1
  Lambdan <- t(Xmat)%*%Xmat + Lambda0t; Lambdaninv = round(ginv(Lambdan), 9)
  B0 = matrix(c(rep(0, 3)), nrow = 1); Lambdans = t(Xmat)%*%Xmat + Lambda0t
  betaTn = ginv(t(Xmat)%*%Xmat + Lambda0t) %*% (t(Xmat)%*% ST[,1:3] + Lambda0t %*% B0)
  betaT = rmatrixnorm(betaTn, Lambdaninv, as.symmetric.matrix(Sig))      # dim nxk (3 x 2), nxn (3x3), kxk (2x2)         
  
  holdalpha1[sim]=betaT[1]; holdbeta0[sim]=betaT[2]; holdbeta1[sim]=betaT[3]
  
  if(Xfit == TRUE){
    
    Xmat = rep(1, n)
    Lambda0t = matrix(c(rep(0, 1)), nrow = 1); diag(Lambda0t) = .1
    theta0S = 0; tauST0 <- holdS[1,1,sim-1]^2
    v<-ginv(Lambda0t+as.numeric(tauST0)*(t(Xmat)%*%Xmat))
    m<-v %*% (tauST0*t(Xmat) %*% ST[,1])  # more generally, diag(t0mat)%*%m0
    betaT<-c(rmvnorm(1,m,v)); holdalpha1[sim]=betaT[1]
    
    Xmat = cbind(1, X) # n by k
    Lambda0t = matrix(c(rep(0, 4)), nrow = 2); diag(Lambda0t) = .1;   Lambdan <- t(Xmat)%*%Xmat + Lambda0t # put zeroes in
    Lambdaninv = round(ginv(Lambdan), 9);  B0 = matrix(c(rep(0, 3)), nrow = 1)
    Lambdans = t(Xmat)%*%Xmat + Lambda0t
    betaTn = ginv(t(Xmat)%*%Xmat + Lambda0t) %*% (t(Xmat)%*% ST[,1:3])
    betaT = rmatrixnorm(betaTn, Lambdaninv, as.symmetric.matrix(Sig))      # dim nxk (3 x 2), nxn (3x3), kxk (2x2)         
    
    holdalpha1[sim]=betaT[1,1];  holdpsi2[sim]=betaT[2,1]
    holdbeta0[sim]=betaT[1,2];   holdomega1[sim]=betaT[2,2]
    holdbeta1[sim]=betaT[1,3];   holdomega2[sim]=betaT[2,3]
  }
  
  tmp1<-(ST[,1])-Xmat*betaT[1]
  tmp2<-(ST[,2])-Xmat*betaT[2]
  tmp3<-(ST[,3])-Xmat*betaT[3]
  
  if(Xfit == TRUE){
    tmp1<-(ST[,1])-Xmat%*%betaT[,1]
    #tmp1<-(ST[,1])-Xmat[,1]*holdalpha1[sim]
    
    tmp2<-(ST[,2])-Xmat%*%betaT[,2]
    tmp3<-(ST[,3])-Xmat%*%betaT[,3]
  }
  
  resid<-cbind((tmp1),(tmp2),(tmp3))
  
  ##s1 
  s1 = rinvgamma(1, shape = a + n/2, scale = (sum(tmp1^2)/2 + b))  
  s2 = rinvgamma(1, shape = a + n/2, scale = (sum(tmp2^2)/2 + b))  
  s3 = rinvgamma(1, shape = a + n/2, scale = (sum(tmp3^2)/2 + b))  
  
  holdS[,,sim] = S = diag(c(sqrt(s1), sqrt(s2), sqrt(s3)))

  # normal priors for coefficients
  holdmu[,sim] = c(holdalpha1[sim], holdbeta0[sim], holdbeta1[sim]);  holdR[,,sim] = R = trueR

  R = holdR[,,sim-1]
  r23 = holdR[2,3,sim-1]; r13 = holdR[1,3,sim]; r12 = holdR[1,2,sim-1]
  
   U13 = (-2*r12*r23 - sqrt(4*r12^2*r23^2 + 4*(1 - r23^2 - r12^2)))/(-2)
   L13 = (-2*r12*r23 + sqrt(4*r12^2*r23^2 + 4*(1 - r23^2 - r12^2)))/(-2)
    if(is.na(L13)|is.na(U13)) next
  
  low13=ceiling(100*max(min(L13,U13))); up13=min(floor(100*max(L13,U13)))
  d13=up13-low13+1
  
  fr13<-matrix(rep(0,d13*4),d13,4)
  for (k in low13:up13){
    r<-  k/100; R[1,3]=r; R[3,1]=r; fr13[(k-low13+1),1]=r
    if(condindfit == T) R[1,2] = R[2,1] = R[2,3] * R[1,3]
    
    summand = apply(resid, 1, function(resid) t(resid)%*%ginv(S %*% R %*% S)%*%(resid) )
    
    j3 = sum(summand)
    
    fr13[(k-low13+1),2]<-fdelt(n, R, j3)
  }
  
  fr13<-fr13[!is.na(fr13[,2]),]
  fr13<-fr13[!is.infinite(fr13[,2]),]
  fr13<-fr13[complete.cases(fr13),]
  fr13[,2]=fr13[,2]-median(fr13[,2])
  fr13[,2]=exp(fr13[,2])
  fr13<-fr13[!is.infinite(fr13[,2]),]
  m<-which(fr13[,2]==max(fr13[,2]))[1]; fr13[m,]
  x1=max(1, (m-20)); x2=min(length(fr13[,1]), (m+20))
  
  low13=1000*fr13[x1,1]; up13=1000*fr13[x2,1]
  d13=up13-low13+1
  
  fr13<-matrix(rep(0,d13*4),d13,4)
  for (k in low13:up13){
    r<-  k/1000; R[1,3]=r; R[3,1]=r; fr13[(k-low13+1),1]=r
    if(condindfit == T) R[1,2] = R[2,1] = R[2,3] * R[1,3]
    
    summand = apply(resid, 1, function(resid) t(resid)%*%ginv(S %*% R %*% S)%*%(resid) )
    
    j3 = sum(summand)
    
    fr13[(k-low13+1),2]<-fdelt(n, R, j3)
  }
  fr13<-fr13[!is.na(fr13[,2]),]
  fr13<-fr13[!is.infinite(fr13[,2]),]
  fr13<-fr13[complete.cases(fr13),]
  fr13[,2]=fr13[,2]-median(fr13[,2])
  fr13[,2]=exp(fr13[,2])
  fr13<-fr13[!is.infinite(fr13[,2]),]
  for (k in 1: length(fr13[,1])){
    fr13[k,3]<-fr13[k,2]/(sum(fr13[,2]))
  }
  
  for (k in 1:length(fr13[,1])){
    fr13[k,4]<-sum(fr13[1:k,3])
  }
  
  u<-runif(1,0,1)
  if (u<fr13[1,4]){
    r13=fr13[1,1]
  }
  if (u>fr13[length(fr13[,1]),4]){
    r13=fr13[length(fr13[,1]),1]
  }
  if (u>=fr13[1,4] & u<=fr13[length(fr13[,1]),4]){
    nr<-nearest(fr13[,4], u)
    r13=fr13[nr,1]
  }
  r13=r13[1]
  R[1,3]=r13; R[3,1]=r13
  
  ### r23 rt
  r23 = R[2,3]; r13 = R[1,3]; r12 = R[1,2]
  U23 = (-2*r12*r13-sqrt(4*r12^2*r13^2 + 4*(1 - r13^2 - r12^2)))/(-2)
  L23 = (-2*r12*r13+sqrt(4*r12^2*r13^2 + 4*(1 - r13^2 - r12^2)))/(-2)
  if(is.na(L23)|is.na(U23)) next
  
  low23=ceiling(100*max(min(L23,U23))); up23=min(floor(100*max(L23,U23)))
  d23=up23-low23+1
  
  fr23<-matrix(rep(0,d23*4),d23,4)
  for (k in low23:up23){
    r<-  k/100; R[2,3]=r; R[3,2]=r; fr23[(k-low23+1),1]=r
    if(condindfit == T) R[1,2] = R[2,1] = R[2,3] * R[1,3]
    summand = apply(resid, 1, function(resid) t(resid)%*%ginv(S %*% R %*% S)%*%(resid) )
    j3 = sum(summand)

    fr23[(k-low23+1),2]<-fdeltBeta(n, R, j3)
 #fr23[(k-low23+1),2]<-fdelt(n, R, j3)
  }
  
  fr23<-fr23[!is.na(fr23[,2]),]
  fr23<-fr23[!is.infinite(fr23[,2]),] 
  fr23<-fr23[complete.cases(fr23),]
  fr23[,2]=fr23[,2]-median(fr23[,2]); fr23[,2]=exp(fr23[,2])
  fr23<-fr23[!is.infinite(fr23[,2]),]; m<-which(fr23[,2]==max(fr23[,2]))[1]; fr23[m,]
  x1=max(1, (m-10)); x2=min(length(fr23[,1]), (m+10))
  
  low23=1000*fr23[x1,1]; up23=1000*fr23[x2,1]
  d23=up23-low23+1
  
  fr23<-matrix(rep(0,d23*4),d23,4)
  for (k in low23:up23){
    r <-  k/1000; R[2,3]=r; R[3,2]=r; fr23[(k-low23+1),1] = r
    if(condindfit == T) R[1,2] = R[2,1] = R[2,3] * R[1,3]
    
    summand = apply(resid, 1, function(resid) t(resid)%*%ginv(S %*% R %*% S)%*%(resid) )
    
    j3 = sum(summand)
    
    fr23[(k-low23+1),2]<-fdeltBeta(n, R, j3)
fr23[(k-low23+1),2]<-fdelt(n, R, j3)
  }
  fr23<-fr23[!is.na(fr23[,2]),]
  fr23<-fr23[!is.infinite(fr23[,2]),] 
  fr23<-fr23[complete.cases(fr23),]
  fr23[,2]=fr23[,2]-median(fr23[,2])
  fr23[,2]=exp(fr23[,2])
  fr23<-fr23[!is.infinite(fr23[,2]),] 
  for (k in 1: length(fr23[,1])){
    fr23[k,3]<-fr23[k,2]/(sum(fr23[,2]))
  }
  
  for (k in 1:length(fr23[,1])){
    fr23[k,4]<-sum(fr23[1:k,3])
  }
  
  u<-runif(1,0,1)
  if (u<fr23[1,4]){
    r23=fr23[1,1]
  }
  if (u>fr23[length(fr23[,1]),4]){
    r23=fr23[length(fr23[,1]),1]
  }
  if (u>=fr23[1,4] & u<=fr23[length(fr23[,1]),4]){
    nr<-nearest(fr23[,4], u)
    r23=fr23[nr,1]
  }
  r23=r23[1]
  R[2,3]=r23; R[3,2]=r23
  
  R[1,2] = R[2,1] = R[1,3] * R[2,3]
  
  if(condindfit == FALSE){
    U12 = (-2*r13*r23-sqrt(4*r13^2*r23^2 + 4*(1 - r23^2 - r13^2)))/(-2)
    L12 = (-2*r13*r23+sqrt(4*r13^2*r23^2 + 4*(1 - r23^2 - r13^2)))/(-2)
    if(is.na(L12)|is.na(U12)) next
    low12=ceiling(100*max(min(L12,U12), -1))
    up12=floor(min(100*max(L12,min(U12),1)))
    d12=up12-low12+1
    
    fr12<-matrix(rep(0,d12*4),d12,4)
    for (k in low12:up12){
      r<- k/100;  R[1,2]=r; R[2,1]=r
      fr12[(k-low12+1),1]=r
      
      summand=NULL
      for(i in 1:n){
        summand[i]=t(resid[i,])%*%ginv(S %*% R %*% S)%*%(resid[i,])
      }
      j3 = sum(summand)
      
      fr12[(k-low12+1),2]<-fdelt(n, R, j3)
      
    }
    fr12=matrix(fr12,ncol=4)
    fr12<-fr12[!is.na(fr12[,2]),] 
    fr12=matrix(fr12,ncol=4)
    fr12<-fr12[!is.infinite(fr12[,2]),] 
    fr12=matrix(fr12,ncol=4)
    fr12<-fr12[complete.cases(fr12),]
    fr12=matrix(fr12,ncol=4)
    fr12[,2]=fr12[,2]-median(fr12[,2])
    fr12=matrix(fr12,ncol=4)
    fr12[,2]=exp(fr12[,2])
    fr12=matrix(fr12,ncol=4)
    fr12<-fr12[!is.infinite(fr12[,2]),]  
    fr12=matrix(fr12,ncol=4)
    m<-which(fr12[,2]==max(fr12[,2]))[1]
    x1=max(1, (m-10)); x2=min(length(fr12[,1]), (m+10))
    low12=1000*fr12[x1,1]; up12=1000*fr12[x2,1]
    d12=up12-low12+1
    
    fr12<-matrix(rep(0,d12*4),d12,4)
    for (k in low12:up12){
      r<- k/1000;  R[1,2]=r;  R[2,1]=r
      fr12[(k-low12+1),1]=r
      
      summand=NULL
      for(i in 1:n){
        summand[i]=t(resid[i,])%*%ginv(S %*% R %*% S)%*%(resid[i,])
      }
      j3 = sum(summand)
      
      fr12[(k-low12+1),2]<-fdelt(n, R, j3)
    }
    fr12=matrix(fr12,ncol=4)
    fr12=matrix(fr12,ncol=4)
    fr12<-fr12[!is.na(fr12[,2]),] 
    fr12=matrix(fr12,ncol=4)
    fr12<-fr12[!is.infinite(fr12[,2]),] 
    fr12=matrix(fr12,ncol=4)
    fr12<-fr12[complete.cases(fr12),]
    fr12=matrix(fr12,ncol=4)
    fr12[,2]=fr12[,2]-median(fr12[,2])
    fr12=matrix(fr12,ncol=4)
    fr12[,2]=exp(fr12[,2])
    fr12=matrix(fr12,ncol=4)
    fr12<-fr12[!is.infinite(fr12[,2]),]  
    fr12=matrix(fr12,ncol=4)
    for (k in 1: length(fr12[,1])){
      fr12[k,3]<-fr12[k,2]/(sum(fr12[,2]))
    }
    
    for (k in 1:length(fr12[,1])){
      fr12[k,4]<-sum(fr12[1:k,3])
    }
    u<-runif(1,0,1)
    if (u<fr12[1,4]){
      r12=fr12[1,1]
    }
    if (u>fr12[length(fr12[,1]),4]){
      r12=fr12[length(fr12[,1]),1]
    }
    if (u>=fr12[1,4] & u<=fr12[length(fr12[,1]),4]){
      nr<-nearest(fr12[,4], u); r12=fr12[nr,1]
    }
    R[1,2] = R[2,1] = r12=r12[1]
  }
  }
  
  if(impute == FALSE){
  
  tauST1<-holdS[3,3,sim-1]
  tauST0<-holdS[2,2,sim-1]
  tauSS0<-holdS[1,1,sim-1]
  XmatS = rep(1, n/2)
  
  v<-ginv(SIG0Smat[1,1]+as.numeric(tauST0)*(t(XmatS)%*%XmatS))
  m<-v%*%(tauST0*t(XmatS)%*%ST[trt == 0,2])  # more generally, diag(t0mat)%*%m0
  betaT<-c(rmvnorm(1,m,v))
  holdbeta0[sim]=betaT[1]
  tmp2<-(ST[trt == 0,2])-XmatS*holdbeta0[sim]
  
  v<-ginv(SIG0Smat[1,1]+as.numeric(tauST1)*(t(XmatS)%*%XmatS))
  m<-v%*%(tauST1*t(XmatS)%*%ST[trt == 1,3])  # more generally, diag(t0mat)%*%m0
  betaT<-c(rmvnorm(1,m,v))
  holdbeta1[sim]=betaT[1]
  tmp3<-(ST[trt == 1,3])-XmatS*holdbeta1[sim]

  v<-ginv(SIG0Smat[1,1]+as.numeric(tauSS0)*(t(XmatS)%*%XmatS))
  m<-v%*%(tauSS0*t(XmatS)%*%ST[trt == 1,1])  # more generally, diag(t0mat)%*%m0
  betaT<-c(rmvnorm(1,m,v))
  holdalpha1[sim]=betaT[1]
  tmp1<-(ST[trt == 1,1])-XmatS*holdalpha1[sim]
  
 if(TRUE){ 
  if(Xfit == FALSE){
  Sigmats = Sig[c(1, 3), c(1, 3)] # 3 x 3
  #Sigmats = holdS[,,1]%*%holdR[,,1]%*%holdS[,,1]
  
  Lambda0t = matrix(c(rep(0, 1)), nrow = 1); diag(Lambda0t) = .1
  Lambdan <- t(XmatS)%*%XmatS + Lambda0t # put zeroes in
  Lambdaninv = round(ginv(Lambdan), 9)
  B0 = matrix(c(rep(0, 2)), nrow = 1)
  Lambdans = t(XmatS)%*%XmatS + Lambda0t
  betaTn = ginv(t(XmatS)%*%XmatS + Lambda0t) %*% (t(XmatS)%*% ST[trt==1, c(1,3)] + Lambda0t %*% B0)
  betaT = rmatrixnorm(betaTn, Lambdaninv, as.symmetric.matrix(Sigmats))      # dim nxk (3 x 2), nxn (3x3), kxk (2x2)         
  
  holdbeta1[sim]=betaT[,2]
  holdalpha1[sim]=betaT[,1] 
  tmp1<-(ST[trt == 1,1])-XmatS*betaT[,1]
  tmp3<-(ST[trt == 1,3])-XmatS*betaT[,2]
  }
   
   if(Xfit == TRUE){
     Xmat = cbind(1, X[trt==0]) # n by k
     
     theta0S = c(0,0)
     v<-ginv(as.numeric(SIG0Smat)+as.numeric(holdS[2,2,sim-1])*(t(Xmat)%*%Xmat))
     m<-v%*%((diag(SIG0Smat)*theta0S)+holdS[2,2,sim-1]*t(Xmat)%*%ST[trt == 0,2])  # more generally, diag(t0mat)%*%m0
     betaT<-c(rmvnorm(1,m,v))
     holdbeta0[sim]=betaT[1]
     holdomega1[sim]=betaT[2]
     
     tmp2<-(ST[trt == 0,2])-Xmat%*%betaT
     
     Xmat = cbind(1, X[trt==1]) # n by k
     
     
     v<-ginv(as.numeric(SIG0Smat)+as.numeric(holdS[3,3,sim-1])*(t(Xmat)%*%Xmat))
     m<-v%*%((diag(SIG0Smat)*theta0S)+holdS[3,3,sim-1]*t(Xmat)%*%ST[trt == 1,3])  # more generally, diag(t0mat)%*%m0
     betaT<-c(rmvnorm(1,m,v))
     holdbeta1[sim]=betaT[1]
     holdomega2[sim]=betaT[2]
     
     tmp3<-(ST[trt == 1,3])-Xmat%*%betaT
     
     v<-ginv(as.numeric(SIG0Smat)+as.numeric(holdS[1,1,sim-1])*(t(Xmat)%*%Xmat))
     m<-v%*%((diag(SIG0Smat)*theta0S)+holdS[1,1,sim-1]*t(Xmat)%*%ST[trt == 1,1])  # more generally, diag(t0mat)%*%m0
     betaT<-c(rmvnorm(1,m,v))
     holdalpha1[sim]=betaT[1]
     holdpsi2[sim]=betaT[2]
     
     tmp1<-(ST[trt == 1,1])-Xmat%*%betaT
   }
}
  
  holdmu[,sim] = c(holdalpha1[sim], holdbeta0[sim], holdbeta1[sim])

  # inv gamma priors for S
  s1 =   rinvgamma(1, shape = a + n/4, scale = (sum(tmp1^2)/2 + b))
  s2 =   rinvgamma(1, shape = a + n/4, scale = (sum(tmp2^2)/2 + b))  
  s3 =   rinvgamma(1, shape = a + n/4, scale = (sum(tmp3^2)/2 + b))  
  
  holdS[,,sim] = S = diag(c(sqrt(s1), sqrt(s2), sqrt(s3)))
  
  ###r13 rho11
  R = holdR[,,sim-1]
  L13 = -.999; U13 = .999
  if(is.na(L13)|is.na(U13)) next
  
  low13=ceiling(100*max(min(L13,U13))); up13=min(floor(100*max(L13,U13)))
  d13=up13-low13+1
  resid = cbind(tmp1, tmp3)
  fr13<-matrix(rep(0,d13*4),d13,4)
  for (k in low13:up13){
    r <-  k/100; R[1,3]=r; R[3,1]=r; fr13[(k-low13+1),1]=r
    if((det(S[c(1,3), c(1,3)] %*% 
               R[c(1,3), c(1,3)] %*% S[c(1,3), c(1,3)])) < 0) next
    
    summand = apply(resid, 1, function(resid) resid%*%ginv(S[c(1,3), c(1,3)] %*% 
         R[c(1,3), c(1,3)] %*% S[c(1,3), c(1,3)])%*%resid)
    
    j3 = sum(summand, na.rm = T)
    fr13[(k-low13+1),2]<-fdelt(n/2, R[c(1,3), c(1,3)], j3)
  }
  
  fr13<-fr13[!is.na(fr13[,2]),]
  fr13<-fr13[!is.infinite(fr13[,2]),] 
  fr13<-fr13[complete.cases(fr13),]
  fr13[,2]=fr13[,2]-median(fr13[,2])
  fr13[,2]=exp(fr13[,2])
  fr13<-fr13[!is.infinite(fr13[,2]),]  
  m<-which(fr13[,2]==max(fr13[,2]))[1]; fr13[m,]
  x1=max(1, (m-10)); x2=min(length(fr13[,1]), (m+10))
  
  low13=1000*fr13[x1,1]; up13=1000*fr13[x2,1]
  d13=up13-low13+1
  if(is.na(d13)) next
  
  fr13<-matrix(rep(0,d13*4),d13,4)
  for (k in low13:up13){
    r <-  k/1000; R[1,3]=r; R[3,1]=r; fr13[(k-low13+1),1] = r
    summand = apply(resid, 1, function(resid) resid%*%ginv(S[c(1,3), c(1,3)] %*% 
                                                             R[c(1,3), c(1,3)] %*% S[c(1,3), c(1,3)])%*%resid)
    
    j3 = sum(summand, na.rm = T)
    
    fr13[(k-low13+1),2]<-fdelt(n/2, R[c(1,3), c(1,3)], j3)
  }
  fr13<-fr13[!is.na(fr13[,2]),]
  fr13<-fr13[!is.infinite(fr13[,2]),] 
  fr13<-fr13[complete.cases(fr13),]
  fr13[,2]=fr13[,2]-median(fr13[,2])
  fr13[,2]=exp(fr13[,2])
  fr13<-fr13[!is.infinite(fr13[,2]),] 
  if(nrow(fr13) == 0) next
  
  for (k in 1: length(fr13[,1])){
    fr13[k,3]<-fr13[k,2]/(sum(fr13[,2]))
  }
  
  for (k in 1:length(fr13[,1])){
    fr13[k,4]<-sum(fr13[1:k,3])
  }
  
  u<-runif(1,0,1)
  if (u<fr13[1,4]){
    r13=fr13[1,1]
  }
  if (u>fr13[length(fr13[,1]),4]){
    r13=fr13[length(fr13[,1]),1]
  }
  if (u>=fr13[1,4] & u<=fr13[length(fr13[,1]),4]){
    nr<-nearest(fr13[,4], u)
    r13=fr13[nr,1]
  }
  
  r13=r13[1]
  R[1,3]=r13; R[3,1]=r13

  ## rt
  R[2,3] = R[3,2] = r23 = runif(1, -1, 1)

  R[2,3] = R[3,2] = r23 = holdR[2,3,sim] = rBeta_ab(1, shape1 = 5, shape2 = 6, a = -0.4, b = 1)
  #R[2,3] = R[3,2] = r23 = holdR[2,3,sim] = 0.22
  #R[2,3] = R[3,2] = r23 = holdR[2,3,sim] = 0.85
  #R[2,3] = holdR[2,3,1]
  
  ###r12 rt
  R[1,2] = R[2,1] = R[1,3] * R[2,3]
  
  if(condindfit == FALSE){
    U12 = (-2*r13*r23-sqrt(4*r13^2*r23^2 + 4*(1 - r23^2 - r13^2)))/(-2)
    L12 = (-2*r13*r23+sqrt(4*r13^2*r23^2 + 4*(1 - r23^2 - r13^2)))/(-2)
    if(is.na(U12) | is.na(L12)) next
    R[1,2] = R[2,1] = runif(1,L12,U12)
  }

  holdR[,,sim] = R
}
  

if(any(eigen(R)$values<0)) next
if(any(abs(R)>1)) next
holdR[,,sim] = R
slope[sim]=(holdR[1,3,sim]*holdS[3,3,sim]-holdR[1,2,sim]*holdS[2,2,sim])/(holdS[1,1,sim])

#slope[sim] = 0.55
int[sim]=(holdmu[3,sim]-holdmu[2,sim])-slope[sim]*(holdmu[1,sim])
holdmu[,sim] = c(holdalpha1[sim], holdbeta0[sim], holdbeta1[sim])

prob_t1t0[sim,] = surrogprob(0, S1, gamma0 = int[sim], gamma1 = slope[sim], 
                  rho10 = R[1,2], rhot = R[2,3], rho11 = R[1,3], sigs1 = holdS[1,1,sim],
                  sigt0 = holdS[2,2,sim], sigt1 = holdS[3,3,sim])
if(Xfit == TRUE){
  holdmu1[,sim] = c(holdalpha1[sim] + holdpsi2[sim], holdbeta0[sim] + holdomega1[sim], holdbeta1[sim] + holdomega2[sim])
  int1[sim]=(holdmu1[3,sim]-holdmu1[2,sim])-slope[sim]*(holdmu1[1,sim])
}

s =  seq(min(samp[,1], -0.01), max(samp[,1]), length.out = n) # grid should stay same every time
intgrid = probgrid = probgrid2 = probgridint = weightgrid = weightgridclosed = Xgrid = sgrid = matrix(s, length(s), length(s))

if(Xfit == TRUE){ 
  l = 1
  eta1h = mean(ST[,1]); eta2h = mean(ST[,2]); eta3h = mean(ST[,3]); eta4h = mean(X) # mean terms
  tau1h = sd(ST[,1]); tau2h = sd(ST[,2]); tau3h = sd(ST[,3]); tau4h = sd(X) # diagonal sd terms
  pi10h = cor(ST[,1], ST[,2]); pi11h = cor(ST[,1], ST[,3]); pi1Xh = cor(ST[,1], X)
  piTh = cor(ST[,2], ST[,3]); piX0h = cor(ST[,2], X)
  piX1h = cor(ST[,3], X)
  
  Xgrid[,1] = X 
  for(l in 2:ncol(sgrid)){
    Xgrid[,l] = c(X[l:(n)], X[1:(l - 1)])
  }
  
  for(l in 1:ncol(sgrid)){
    
    Xpermute = Xgrid[,l]
    intgrid[,l] = intx = ((holdbeta1[sim] + holdomega2[sim]*Xpermute) -  (holdbeta0[sim] + holdomega1[sim]*Xpermute))-slope[sim] * (
                              holdalpha1[sim] + holdpsi2[sim]*Xpermute )
    sumSX = 0*Xpermute; xindex = 1
    
    for(xindex in 1:nrow(sgrid)){                       
      sumSX[xindex] = exp(-1/(2 * S[1,1]^2) * (sgrid[xindex, l] - (c(1, Xgrid[xindex, l]) %*% c(holdalpha1[sim], holdpsi2[sim]) ) )^2) # gives correct answer
    } # gives same as closed form
    
    
    weightgrid[,l] = sumSX  # /sum(sumSX)
    #weightgridclosed[,l] = dnorm(Xgrid[,l], mean = eta4 + tau1^(-1)*tau4*pi1X*(s - eta1), sd = sqrt(tau4^2 - tau4^2*pi1X^2)) # X|S dist
    weightgridclosed[,l] = dnorm(Xgrid[,l], mean = eta4h + tau1h^(-1)*tau4h*pi1Xh*(s - eta1h), sd = sqrt(tau4h^2 - tau4h^2*pi1Xh^2)) # X|S dist
  }
  
  weightgrid = weightgrid / rowSums(weightgrid) # weights are normalized so that each weight for a given s adds to 1
  weightgridclosed = weightgridclosed / rowSums(weightgridclosed) # weights are normalized so that each weight for a given s addss to 1
  
  xindex = 1; l = 1
  for(l in 1:nrow(sgrid)){
    
    for(xindex in 1:nrow(sgrid)){
      probgrid[xindex, l] = ( intgrid[xindex, l] + slope[sim]*sgrid[xindex, l]) * weightgrid[xindex, l] #* dnorm(Xpermute[xindex], mean = 0.5, sd = 1) # times f(x)
    }
  }

  for(l in 1:nrow(sgrid)){
    for(xindex in 1:nrow(sgrid)){
      probgrid2[xindex, l] = ( intgrid[xindex, l] + slope[sim]*sgrid[xindex, l]) * weightgridclosed[xindex, l] #* dnorm(Xpermute[xindex], mean = 0.5, sd = 1) # times f(x)
    }
  }
  if(sim <= 5){
    plot(s, rowSums(probgrid)); lines(s, int[sim] + slope[sim]*s, col = "red"); if(Xfit == TRUE){lines(s, int1[sim] + slope[sim]*s, col = "blue")}
    plot(s, rowSums(probgrid2)); lines(s, int[sim] + slope[sim]*s, col = "red"); if(Xfit == TRUE){lines(s, int1[sim] + slope[sim]*s, col = "blue")}
  }
  
  intmarginalize[sim] = lm((rowSums(probgrid)) ~ s)$coef[1]
  slopemarginalize[sim] = lm((rowSums(probgrid)) ~ s)$coef[2]
  intmarginalizewtcl[sim] = lm((rowSums(probgrid2)) ~ s)$coef[1]
  slopemarginalizewtcl[sim] = lm((rowSums(probgrid2)) ~ s)$coef[2]
  
}

lm(ST[,3] - ST[,2] ~ ST[,1])

if(Xfit == TRUE){
  A = holdbeta1[sim] - holdbeta0[sim] - slope[sim]*holdalpha1[sim]
  B = holdomega2[sim] - holdomega1[sim] - slope[sim]*holdpsi2[sim]
  eta1h = mean(ST[,1]); eta2h = mean(ST[,2]); eta3h = mean(ST[,3]); eta4h = mean(X) # mean terms
  tau1h = sd(ST[,1]); tau2h = sd(ST[,2]); tau3h = sd(ST[,3]); tau4h = sd(X) # diagonal sd terms
  pi10h = cor(ST[,1], ST[,2]); pi11h = cor(ST[,1], ST[,3]); pi1Xh = cor(ST[,1], X)
  piTh = cor(ST[,2], ST[,3]); piX0h = cor(ST[,2], X)
  piX1h = cor(ST[,3], X)
  
  gamma0star = A + B*(eta4h - tau1h^(-1)*tau4h*pi1Xh*eta1h) # aligns with previous belief
  gamma1star = B*(tau1h^(-1)*tau4h*pi1Xh) + slope[sim] # works with no X effect
  if(sim <= 5){
    plot(s,  gamma0star + gamma1star*s) # corresponds to what I solved for before
    lines(s, int[sim] + slope[sim]*s, col = "red");
    if(Xfit == TRUE) lines(s, int1[sim] + slope[sim]*s, col = "blue");
  }
  intmarginalizeclosed[sim] = gamma0star
  slopemarginalizeclosed[sim] = gamma1star
  
}

print(R)
print(slope[sim])
print(slopemarginalizewtcl[sim])

sim=sim+1
print(sim)

}    

gam<-array(0,c((sim-1),1))
alph<-array(0,c((sim-1),1))

b1<-array(0,c((sim-1),1))
b2<-array(0,c((sim-1),1))

slope[1]=(holdR[1,3,1]*holdS[3,3,1]-holdR[1,2,1]*holdS[2,2,1])/(holdS[1,1,1])
int[1]=(holdmu[3,1]-holdmu[2,1])-slope[1]*(holdmu[1,1])
int1[1]=(holdmu1[3,1]-holdmu1[2,1])-slope[1]*(holdmu1[1,1])

prob_t1t0[1,] = surrogprob(0, S1, gamma0 = int[1], gamma1 = slope[1], 
                           rho10 = holdR[1,2,1], rhot = holdR[2,3,1], rho11 = holdR[1,3,1], sigs1 = holdS[1,1,sim],
                           sigt0 = holdS[2,2,1], sigt1 = holdS[3,3,1])

hist(holdR[1,3,])
hist(holdR[1,2,])
hist(holdR[2,3,])

prob_t1t0l = apply(X = prob_t1t0, MARGIN = 2, FUN = quantile, probs = 0.025, na.rm=T)
prob_t1t0u = apply(X = prob_t1t0, MARGIN = 2, FUN = quantile, probs = 0.975, na.rm=T)

plot(S1, colMeans(prob_t1t0, na.rm = T), xlab = "S1 Value", ylab = "Probability(T1-T0>0|S1)", main = 
       "Treatment Efficacy by Surrogate,\nImputation", ylim = c(0, 1)); lines(S1, prob_t1t0l); lines(S1, prob_t1t0u)

mu<-array(0,c((sim-1),1))
params<-data.frame(beta0SE=numeric(1), beta1SE=numeric(1), alpha0SE=numeric(1),  alpha1SE=numeric(1),
                   psi1SE=numeric(1), omega1SE=numeric(1), psi2SE=numeric(1), omega2SE=numeric(1),
                   sigs0=numeric(1), SEsigs0=numeric(1), sigs1=numeric(1), SEsigs1=numeric(1),
                   sigmat0=numeric(1), SEsigmat0=numeric(1),sigmat1=numeric(1), SEsigmat1=numeric(1), 
                   ps=numeric(1), SEps=numeric(1),p00=numeric(1), SEp00=numeric(1),
                   p01=numeric(1),SEp01=numeric(1),p10=numeric(1),SEp10=numeric(1),p11=numeric(1),
                   SEp11=numeric(1),pt=numeric(1),SEpt=numeric(1), S0l=numeric(1), S0u=numeric(1),
                   S1l=numeric(1), S1u=numeric(1), T0l=numeric(1),T0u=numeric(1), T1l=numeric(1),
                   T1u=numeric(1), psl=numeric(1), PSu=numeric(1),p00l=numeric(1), p00u=numeric(1),
                   p01l=numeric(1),p01u=numeric(1),p10l=numeric(1),p10u=numeric(1),p11l=numeric(1),
                   p11u=numeric(1),ptl=numeric(1),ptu=numeric(1))

params[1]<-sqrt(var(holdbeta0[burnin:(sim-1)]))
params[2]<-sqrt(var(holdbeta1[burnin:(sim-1)]))
#params[3]<-sqrt(var(holdalpha0[burnin:(sim-1)]))
params[4]<-sqrt(var(holdalpha1[burnin:(sim-1)]))
params[5]<-sqrt(var(holdpsi1[burnin:(sim-1)]))
params[6]<-sqrt(var(holdomega1[burnin:(sim-1)]))
params[7]<-sqrt(var(holdpsi2[burnin:(sim-1)]))
params[8]<-sqrt(var(holdomega2[burnin:(sim-1)]))

params[9]<-mean(holdS[1,1,burnin:(sim-1)])
params[10]<-sqrt(var(holdS[1,1,burnin:(sim-1)]))
params[11]<-mean(holdS[2,2,burnin:(sim-1)])
params[12]<-sqrt(var(holdS[2,2,burnin:(sim-1)]))
params[13]<-mean(holdS[3,3,burnin:(sim-1)])
params[14]<-sqrt(var(holdS[3,3,burnin:(sim-1)]))

params[17]<-mean(holdR[1,2,burnin:(sim-1)])
params[18]<-sqrt(var(holdR[1,2,burnin:(sim-1)]))
params[19]<-mean(holdR[1,3,burnin:(sim-1)])
params[20]<-sqrt(var(holdR[1,3,burnin:(sim-1)]))
params[23]<-mean(holdR[2,3,burnin:(sim-1)])
params[24]<-sqrt(var(holdR[2,3,burnin:(sim-1)]))

params[29]<-quantile(holdS[1,1,burnin:(sim-1)],  probs = 0.025,na.rm=TRUE)
params[30]<-quantile(holdS[1,1,burnin:(sim-1)],  probs = 0.975,na.rm=TRUE)
params[31]<-quantile(holdS[2,2,burnin:(sim-1)],  probs = 0.025,na.rm=TRUE)
params[32]<-quantile(holdS[2,2,burnin:(sim-1)],  probs = 0.975,na.rm=TRUE)
params[33]<-quantile(holdS[3,3,burnin:(sim-1)],  probs = 0.025,na.rm=TRUE)
params[34]<-quantile(holdS[3,3,burnin:(sim-1)],  probs = 0.975,na.rm=TRUE)

params[37]<-quantile(holdR[1,2,burnin:(sim-1)],  probs = 0.025,na.rm=TRUE)
params[38]<-quantile(holdR[1,2,burnin:(sim-1)],  probs = 0.975,na.rm=TRUE)
params[39]<-quantile(holdR[1,3,burnin:(sim-1)],  probs = 0.025,na.rm=TRUE)
params[40]<-quantile(holdR[1,3,burnin:(sim-1)],  probs = 0.975,na.rm=TRUE)
params[43]<-quantile(holdR[2,3,burnin:(sim-1)],  probs = 0.025,na.rm=TRUE)
params[44]<-quantile(holdR[2,3,burnin:(sim-1)],  probs = 0.975,na.rm=TRUE)

prentice<-data.frame(bet_s=numeric(1), SEbs=numeric(1),bet_trtt=numeric(1), SEbtt=numeric(1), 
                     bet_trts=numeric(1), SEbts=numeric(1), gam_trt=numeric(1), SEgt=numeric(1),
                     gam_s=numeric(1), SEgs=numeric(1))

pren.post<-data.frame(bet_sp=numeric(1), SEbet_sp=numeric(1),bet_trttp=numeric(1), SEbet_trttp=numeric(1),  
                      bet_trtsp=numeric(1), SEbet_trtsp=numeric(1), gam_trtp=numeric(1), SEgam_trtp=numeric(1), 
                      gam_sp=numeric(1),SEgam_sp=numeric(1))

PS<-data.frame(dat_int=numeric(1), dat_intSE=numeric(1),dat_sl=numeric(1), dat_slSE=numeric(1), 
               mean_int=numeric(1), SEmean_int=numeric(1),L_int=numeric(1), U_int=numeric(1),mean_sl=numeric(1)
               ,SEmean_sl=numeric(1),L_sl=numeric(1),U_sl=numeric(1),
               postdat_int=numeric(1), postdat_intSE=numeric(1), postdat_sl=numeric(1), postdat_slSE=numeric(1),
               mean_int=numeric(1), SEmean_int=numeric(1),L_int=numeric(1), U_int=numeric(1), covslint=numeric(1),covslint1=numeric(1),
               int_coverage=numeric(1), slope_coverage=numeric(1), int1_coverage=numeric(1),
               intmarg_mean=numeric(1), intmarg_SE=numeric(1), intmargl=numeric(1), intmargu=numeric(1), intmarg_coverage=numeric(1),
               intaccept=numeric(1), int1accept=numeric(1), intmargaccept=numeric(1), slopeaccept=numeric(1),
               slopemarg_mean=numeric(1), slopemarg_SE=numeric(1), slopemargl=numeric(1), slopemargu=numeric(1),
               intmarg_meancl=numeric(1), intmarg_SEcl=numeric(1), intmarglcl=numeric(1), intmargucl=numeric(1),
               slopemarg_meancd=numeric(1), slopemarg_SEcd=numeric(1), slopemarglcd=numeric(1), slopemargucd=numeric(1),
               intmarg_meanwtcl=numeric(1), intmarg_SEwtcl=numeric(1), intmarglwtcl=numeric(1), intmarguwtcl=numeric(1),
               slopemarg_meanwtcl=numeric(1), slopemarg_SEwtcl=numeric(1), slopemarglwtcl=numeric(1), slopemarguwtcl=numeric(1)
)

covs<-data.frame(psl=numeric(1), psu=numeric(1),p00l=numeric(1), p00u=numeric(1),
                 p01l=numeric(1),p01u=numeric(1),p10l=numeric(1),p10u=numeric(1),p11l=numeric(1),
                 p11u=numeric(1),ptl=numeric(1),ptu=numeric(1),s0l=numeric(1),s0u=numeric(1)
                 ,s1l=numeric(1),s1u=numeric(1),t0l=numeric(1),t0u=numeric(1),t1l=numeric(1),t1u=numeric(1),
                 psind=numeric(1),p00ind=numeric(1),p01ind=numeric(1),p10ind=numeric(1),p11ind=numeric(1),ptind=numeric(1),
                 s0ind=numeric(1),s1ind=numeric(1),t0ind=numeric(1),t1ind=numeric(1))

plot(holdR[1,2,])
plot(holdR[2,3,])
plot(holdR[1,3,])

covs[1]<-quantile(holdR[1,2,burnin:sim-1],  probs = 0.025, na.rm=T)
covs[2]<-quantile(holdR[1,2,burnin:sim-1],  probs = 0.975, na.rm=T)
covs[3]<-quantile(holdR[1,3,burnin:sim-1],  probs = 0.025, na.rm=T)
covs[4]<-quantile(holdR[1,3,burnin:sim-1],  probs = 0.975, na.rm=T)
covs[7]<-quantile(holdR[2,3,burnin:sim-1],  probs = 0.025, na.rm=T)
covs[8]<-quantile(holdR[2,3,burnin:sim-1],  probs = 0.975, na.rm=T)

covs[13]<-quantile(holdS[1,1,burnin:sim-1],  probs = 0.025, na.rm=T)
covs[14]<-quantile(holdS[1,1,burnin:sim-1],  probs = 0.975, na.rm=T)
covs[15]<-quantile(holdS[2,2,burnin:sim-1],  probs = 0.025, na.rm=T)
covs[16]<-quantile(holdS[2,2,burnin:sim-1],  probs = 0.975, na.rm=T)
covs[17]<-quantile(holdS[3,3,burnin:sim-1],  probs = 0.025, na.rm=T)
covs[18]<-quantile(holdS[3,3,burnin:sim-1],  probs = 0.975, na.rm=T)

covs[21]<-as.numeric((holdR[1,2,1]>covs[1])&(holdR[1,2,1]<covs[2]))
covs[22]<-as.numeric((holdR[1,3,1]>covs[3])&(holdR[1,2,1]<covs[4]))
covs[24]<-as.numeric((holdR[2,3,1]>covs[7])&(holdR[2,3,1]<covs[8]))
covs[27]<-as.numeric((holdS[1,1,1]>covs[13])&(holdS[1,1,1]<covs[14]))
covs[28]<-as.numeric((holdS[2,2,1]>covs[15])&(holdS[2,2,1]<covs[16]))

PS[5]<-mean(int[burnin:sim-1], na.rm=T)
PS[6]<-sqrt(var(int[burnin:sim-1], na.rm=T))
PS[7]<-quantile(int[burnin:sim-1],  probs = 0.025, na.rm=T)
PS[8]<-quantile(int[burnin:sim-1],  probs = 0.975, na.rm=T)
PS[9]<-mean(slope[burnin:sim-1], na.rm=T)
PS[10]<-sqrt(var(slope[burnin:sim-1], na.rm=T))
PS[11]<-quantile(slope[burnin:sim-1],  probs = 0.025, na.rm=T)
PS[12]<-quantile(slope[burnin:sim-1],  probs = 0.975, na.rm=T)

PS[17]<-mean(int1[burnin:sim-1], na.rm=T)
PS[18]<-sqrt(var(int1[burnin:sim-1], na.rm=T))
PS[19]<-quantile(int1[burnin:sim-1],  probs = 0.025, na.rm=T)
PS[20]<-quantile(int1[burnin:sim-1],  probs = 0.975, na.rm=T)
PS[21]<-cov(int[burnin:sim-1],slope[burnin:sim-1])
PS[22]<-cov(int1[burnin:sim-1],slope[burnin:sim-1])
PS[22]<-cov(int1[burnin:sim-1],slope[burnin:sim-1])
PS[23]=as.numeric(PS[7]<int[1]&PS[8]>int[1])
PS[24]=as.numeric(PS[11]<slope[1]&PS[12]>slope[1])
PS[25]=as.numeric(PS[19]<int1[1]&PS[20]>int1[1])
PS[26]<-mean(intmarginalize[burnin:sim-1], na.rm=T)
PS[27]<-sqrt(var(intmarginalize[burnin:sim-1], na.rm=T))
PS[28]<-quantile(intmarginalize[burnin:sim-1],  probs = 0.025, na.rm=T)
PS[29]<-quantile(intmarginalize[burnin:sim-1],  probs = 0.975, na.rm=T)
PS[30]=as.numeric(PS[28]< 0 & PS[29]> 0) # may have to fix
PS[31]=as.numeric(PS[7]<0 & PS[8]>=0)
PS[32]=as.numeric(PS[19]<0 & PS[20]>=0)
PS[33]=as.numeric(PS[28]<0 & PS[29]>=0)
PS[34]=as.numeric(PS[11]>0)
PS[35]<-mean(slopemarginalize[burnin:sim-1], na.rm=T)
PS[36]<-sqrt(var(slopemarginalize[burnin:sim-1], na.rm=T))
PS[37]<-quantile(slopemarginalize[burnin:sim-1],  probs = 0.025, na.rm=T)
PS[38]<-quantile(slopemarginalize[burnin:sim-1],  probs = 0.975, na.rm=T)
PS[39]<-mean(intmarginalizeclosed[burnin:sim-1], na.rm=T)
PS[40]<-sqrt(var(intmarginalizeclosed[burnin:sim-1], na.rm=T))
PS[41]<-quantile(intmarginalizeclosed[burnin:sim-1],  probs = 0.025, na.rm=T)
PS[42]<-quantile(intmarginalizeclosed[burnin:sim-1],  probs = 0.975, na.rm=T)
PS[43]<-mean(slopemarginalizeclosed[burnin:sim-1], na.rm=T)
PS[44]<-sqrt(var(slopemarginalizeclosed[burnin:sim-1], na.rm=T))
PS[45]<-quantile(slopemarginalizeclosed[burnin:sim-1],  probs = 0.025, na.rm=T)
PS[46]<-quantile(slopemarginalizeclosed[burnin:sim-1],  probs = 0.975, na.rm=T)

PS[47]<-mean(intmarginalizewtcl[burnin:sim-1], na.rm=T)
PS[48]<-sqrt(var(intmarginalizewtcl[burnin:sim-1], na.rm=T))
PS[49]<-quantile(intmarginalizewtcl[burnin:sim-1],  probs = 0.025, na.rm=T)
PS[50]<-quantile(intmarginalizewtcl[burnin:sim-1],  probs = 0.975, na.rm=T)
PS[51]<-mean(slopemarginalizewtcl[burnin:sim-1], na.rm=T)
PS[52]<-sqrt(var(slopemarginalizewtcl[burnin:sim-1], na.rm=T))
PS[53]<-quantile(slopemarginalizewtcl[burnin:sim-1],  probs = 0.025, na.rm=T)
PS[54]<-quantile(slopemarginalizewtcl[burnin:sim-1],  probs = 0.975, na.rm=T)

estcoef<-data.frame(beta0mean=numeric(1), beta1mean=numeric(1),omega1mean=numeric(1), alpha0mean=numeric(1),
                    alpha1mean=numeric(1),psi1mean=numeric(1),psi2mean=numeric(1),omega2mean=numeric(1),
                    beta0l=numeric(1),beta0u=numeric(1), beta1l=numeric(1),beta1u=numeric(1),omega1l=numeric(1),
                    omega1u=numeric(1), alpha0l=numeric(1),alpha0u=numeric(1),
                    alpha1l=numeric(1),alpha1u=numeric(1),psi1l=numeric(1),psi1u=numeric(1),
                    psi2l=numeric(1),psi2u=numeric(1),omega2l=numeric(1),omega2u=numeric(1), 
                    beta0ind=numeric(1), beta1ind=numeric(1),omega1ind=numeric(1), alpha0ind=numeric(1),
                    alpha1ind=numeric(1),psi1ind=numeric(1),psi2ind=numeric(1),omega2ind=numeric(1))

estcoef[1]=mean(holdbeta0[burnin:sim-1], na.rm=T)
estcoef[2]=mean(holdbeta1[burnin:sim-1], na.rm=T)
estcoef[3]=mean(holdomega1[burnin:sim-1], na.rm=T)
#estcoef[4]=alpha0mean=mean(holdalpha0[burnin:sim-1], na.rm=T)
estcoef[5]=mean(holdalpha1[burnin:sim-1], na.rm=T)
estcoef[6]=mean(holdpsi1[burnin:sim-1], na.rm=T)
estcoef[7]=mean(holdpsi2[burnin:sim-1], na.rm=T)
estcoef[8]=omega2mean=mean(holdomega2[burnin:sim-1], na.rm=T)
estcoef[9]<-quantile(holdbeta0[burnin:sim-1],  probs = 0.025, na.rm=T)
estcoef[10]<-quantile(holdbeta0[burnin:sim-1],  probs = 0.975, na.rm=T)
estcoef[11]<-quantile(holdbeta1[burnin:sim-1],  probs = 0.025, na.rm=T)
estcoef[12]<-quantile(holdbeta1[burnin:sim-1],  probs = 0.975, na.rm=T)
estcoef[13]<-quantile(holdomega1[burnin:sim-1],  probs = 0.025, na.rm=T)
estcoef[14]<-quantile(holdomega1[burnin:sim-1],  probs = 0.975, na.rm=T)
#estcoef[15]<-quantile(holdalpha0[burnin:sim-1],  probs = 0.025, na.rm=T)
#estcoef[16]<-quantile(holdalpha0[burnin:sim-1],  probs = 0.975, na.rm=T)
estcoef[17]<-quantile(holdalpha1[burnin:sim-1],  probs = 0.025, na.rm=T)
estcoef[18]<-quantile(holdalpha1[burnin:sim-1],  probs = 0.975, na.rm=T)
estcoef[19]<-quantile(holdpsi1[burnin:sim-1],  probs = 0.025, na.rm=T)
estcoef[20]<-quantile(holdpsi1[burnin:sim-1],  probs = 0.975, na.rm=T)
estcoef[21]<-quantile(holdpsi2[burnin:sim-1],  probs = 0.025, na.rm=T)
estcoef[22]<-quantile(holdpsi2[burnin:sim-1],  probs = 0.975, na.rm=T)
estcoef[23]<-quantile(holdomega2[burnin:sim-1],  probs = 0.025, na.rm=T)
estcoef[24]<-quantile(holdomega2[burnin:sim-1],  probs = 0.975, na.rm=T)
estcoef[25]=as.numeric(estcoef[9]<holdmu[2,1]&estcoef[10]>holdmu[2,1])
estcoef[26]=as.numeric(estcoef[11]<holdmu[3,1]&estcoef[12]>holdmu[3,1])
estcoef[27]=as.numeric(estcoef[13]< holdomega1[1] &estcoef[14]> holdomega1[1])
estcoef[28]=as.numeric(estcoef[15]<holdmu[1,1]&estcoef[16]>holdmu[1,1])
estcoef[29]=as.numeric(estcoef[17]<holdmu[1,1]&estcoef[18]>holdmu[1,1])
#estcoef[30]=as.numeric(estcoef[19]<1&estcoef[20]>1)
estcoef[31]=as.numeric(estcoef[21]< holdpsi2[1] &estcoef[22]> holdpsi2[1])
estcoef[32]=as.numeric(estcoef[23]< holdomega2[1]&estcoef[24]>holdomega2[1])


fname <- paste('params','coefsimVaccinemeans.n',n, method, array_id,'.txt',sep="")
write.table(params, file=fname, sep="\t", row.names=F, col.names=T)
fname2 <- paste('PS','coefsimVaccinemeans.n',n, method, array_id,'.txt',sep="")
write.table(PS, file=fname2, sep="\t", row.names=F, col.names=T)
fname3 <- paste('prentice','coefsim4Vaccine.n', n, method, array_id,'.txt',sep="")
###write.table(prentice, file=fname3, sep="\t", row.names=F, col.names=T)
fname4 <- paste('postpren','coefsimVaccinemeans.n',n, array_id,'.txt',sep="")
###write.table(pren.post, file=fname4, sep="\t", row.names=F, col.names=T)
fname5 <- paste('naivemodels','coefsimVaccinemeans.n',n,array_id,'.txt',sep="")
###write.table(naiveresults, file=fname5, sep="\t", row.names=F, col.names=T)
fname6 <- paste('estimatedcoef','coefsimVaccinemeans.n', n, method, array_id,'.txt',sep="")
write.table(estcoef, file=fname6, sep="\t", row.names=F, col.names=T)
fname7 <- paste('covs','coefSimVaccinemeans.n', n, method, array_id,'.txt',sep="")
write.table(covs, file=fname7, sep="\t",  col.names=T)

plot(int)
plot(int1)
plot(slope)
plot(intmarginalize)
plot(slopemarginalize)

