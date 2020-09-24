library(lme4); library(MASS); library(MASS); library(mvtnorm); library(MCMCpack); library(LaplacesDemon); library(GenKern)
library(ppcor)

array_id = 1; iter = 1
set.seed(12)
SIM = 300
burnin = 2
a = b = 0.1
condindfit = T
difference = F
impute = F
Xfit = TRUE
n = 100
#allsamp = mvrnorm(n, mu = c(delta1, delta2, delta3, delta4, delta5), R) # 5: S1, T0, T1, X, A
#allsamp = allsamp + cbind(psi2 * X, omega1 * X, omega2 * X)

phi1 = 4

phi4 = 0.1
phi7 = 3

phi2 = 0
phi3 = 0
phi5 = 0.4

phi6 = 0.2
phi8 = 0.4
phi9 = 0.3

mu = c(phi1, phi4, phi4) # same intercept

tau = 0.05 # represents ei
allsamp1 = mvrnorm(n, mu = mu, diag(c(tau, tau, tau)))
allsamp2 = mvrnorm(n, mu = mu, diag(c(tau, tau, tau)))
allsamp3 = mvrnorm(n, mu = mu, diag(c(tau, tau, tau)))
allsamp4 = mvrnorm(n, mu = mu, diag(c(tau, tau, tau)))

# CI of RE conditional on X, A scale
psi11 = 0.45
psiT = 0.05
psi10 = 0.09

eta4 = 24
eta5 = 5
pixa = 0.55

R<-matrix(rep(0.22,3*3),3,3); R[1,2] = R[2,1] = psi10;  R[1,3] = R[3,1] <-psi11; 
R[2,3] = R[3,2] <- psiT
R[1,1] = 1
eigen(R)
uitrue = mvrnorm(n, c(mu[1], 0, 0), R)

# S1 T01 T11 T02 T12 T03 T13 T04 T14
Y = cbind(uitrue[,1], allsamp1[,2:3] + uitrue[,2:3], allsamp2[,2:3] + uitrue[,2:3],
          allsamp3[,2:3] + uitrue[,2:3], allsamp4[,2:3] + uitrue[,2:3])
colMeans(Y)

XA = mvrnorm(n, mu = c(eta4, eta5), matrix(c( sqrt(2), 2*sqrt(0.7)*pixa,
                                              2*sqrt(0.7)*pixa, 0.65), nrow = 2) )
A = XA[,2]; X = XA[,1]
coef(lm(Y[,9] - Y[,8] ~  A + I(A^2) + Y[,1] + X))

# add in effect of time and baseline
Y[,2:3] = X
beta_t = 0; beta_t1 = 3.5
Y[,4:5] = Y[,4:5] + beta_t; Y[,6:7] = Y[,6:7] + 2*beta_t; Y[,8:9] = Y[,8:9] + 3*beta_t
# interaction of time and treatment
Y[,5] = Y[,5] + beta_t1; Y[,7] = Y[,7] + 2*beta_t1; Y[,9] = Y[,9] + 3*beta_t1
colMeans(Y)
var(Y[,4:5])

# add in effect of age afterward
phi6 = beta_a = 0.45
phi9 = -1.0
beta_a2 = -0.2
Y[, 4:9] = Y[, 4:9] + beta_a * A + beta_a2 * A^2
Y[, c(5, 7, 9)] = Y[, c(5, 7, 9)] + (phi9 - phi6) * A # interaction with treatment

# subgroup effect, treatment works less well over age of 7
#Y[Age > 7, 4:9] = Y[Age > 7, c(5, 7, 9)] - 3.5

var(Y[,4:5])
colMeans(Y)

cor(Y[,4], Y[,5]) - cor(Y[,1], Y[,4])/cor(Y[,1], Y[,5])
slope = (sd(Y[,5])*cor(Y[,1], Y[,5]) - sd(Y[,4])*cor(Y[,1], Y[,4]))/
  sd(Y[,1])

# add in the effect of X
phi5 = 1.14
phi8 = 1.14
Y[, c(4, 6, 8)] = Y[, c(4, 6, 8)] + phi5 * X
Y[, c(5, 7, 9)] = Y[, c(5, 7, 9)] + phi8 * X

colMeans(Y)

ST = Y[, c(1, 6:7)]
ST = Y[, c(1, 8:9)]

var(ST)
cov(ST, X)
colMeans(ST)
# very negative covariance of ST and X -> taking the difference increases variance of ST
if(difference == TRUE) ST[,2:3] = ST[,2:3] - X

cor(ST[,2], ST[,3]) - cor(ST[,1], ST[,2])/cor(ST[,1], ST[,3])

slope = (sd(ST[,3])*cor(ST[,1], ST[,3]) - sd(ST[,2])*cor(ST[,1], ST[,2]))/sd(ST[,1])
mean(ST[,3]) - mean(ST[,2]) - slope[1]*mean(ST[,1])

cor(ST)
colMeans(ST)

lm(ST ~ X)

allsamp = ST

##### construct dataset
alltrt<-c(rep(0,iter*n)); for (t in 1:iter){  alltrt[(t*n-(n/2-1)):(n*t)]<-1}

# treatment effect at population/CF level
round(mean(ST[,3] - ST[,2]), 3)
mean(ST[,3] - X) - mean(ST[,2] - X)
lm(ST[,3] - ST[,2] ~ A)
lm(ST[,3] - ST[,2] ~ A + X)
lm(ST[,3] - ST[,2] ~ X)

# interactions for treatment effect and surrogacy in subgroups
sd(ST[,3])
sd(ST[,2])

lm(ST[,3] - ST[,2] ~ ST[,1])
lm(ST[,3] - ST[,2] ~ A + ST[,1])
lm(ST[,3] - ST[,2] ~ A + I(A^2) + ST[,1])
lm(ST[,3] - ST[,2] ~ X + A + I(A^2) + ST[,1])

lm(ST~ X + A + I(A^2))
lm(ST[trt==1,]~ X[trt==1] + A[trt==1] + I(A[trt==1]^2))

Xcent = (X - mean(X))
Acent = (A - mean(A))

lm(ST[,3] - ST[,2] ~ Xcent + Acent + I(Acent^2) + ST[,1])

y = coef(lm(ST[,3] - ST[,2] ~  A + I(A^2) + ST[,1]))
y = coef(lm(ST[,3] - ST[,2] ~  A + ST[,1]))

y[1] + y[2] * 5
y[1] + y[2] * 4
y[1] + y[2] * 6

x = coef(lm(ST[,3] - ST[,2] ~  A + I(A^2) + ST[,1] + X))
x[1] + x[2] * 5 + x[3] * 5^2 + x[5]*12
x[1] + x[2] * 5 + x[3] * 5^2 + x[5]*20
x[1] + x[2] * 5 + x[3] * 5^2 + x[5]*26

x[1] + x[2] * 6 + x[3] * 6^2 + x[5]*20
x[1] + x[2] * 6 + x[3] * 6^2 + x[5]*16

x[1] + x[2] * 2 + x[3] * 2^2 + x[5]*12
x[1] + x[2] * 2 + x[3] * 2^2 + x[5]*20

# partial correlations
# is the correlation between the residuals
# are MVN, the partial correlation equals the conditional correlation.

pcor(cbind(ST[,3] - ST[,2], A, X))$estimate
pcor(cbind(ST[,3] - ST[,2], A, X, ST[,1]))$estimate # cor T1-T0, S1 | A, X
pcor(cbind(ST[,3] - ST[,2], A, ST[,1]))$estimate
pcor(cbind(ST[,3] - ST[,2], X, ST[,1]))$estimate

condaxt1 = pcor(cbind(ST[,3], A, X, ST[,1], A^2))$estimate # cor T1, S1 | A, X
condaxt0 = pcor(cbind(ST[,2], A, X, ST[,1], A^2))$estimate # cor T0, S1 | A, X
condaxt = pcor(cbind(ST[,2], A, X, ST[,3], A^2))$estimate # cor T0, T1 | A, X

condaxt[1,4] - condaxt0[1,4]/condaxt1[1,4] # CI conditional on A, X
condaxt1[1,4] - condaxt0[1,4]

condaxt1 = pcor(cbind(ST[,3], A, X, ST[,1]))$estimate # cor T1, S1 | A, X
condaxt0 = pcor(cbind(ST[,2], A, X, ST[,1]))$estimate # cor T0, S1 | A, X
condaxt = pcor(cbind(ST[,2], A, X, ST[,3]))$estimate # cor T0, T1 | A, X

condaxt[1,4] - condaxt0[1,4]/condaxt1[1,4] # CI conditional on A, X

condat1 = pcor(cbind(ST[,3], A, ST[,1]))$estimate # cor T1, S1 | A
condat0 = pcor(cbind(ST[,2], A, ST[,1]))$estimate # cor T0, S1 | A
condat = pcor(cbind(ST[,2], A, ST[,3]))$estimate # cor T0, T1 | A

condat1 = pcor(cbind(ST[,3], A, ST[,1], A^2))$estimate # cor T1, S1 | A
condat0 = pcor(cbind(ST[,2], A, ST[,1], A^2))$estimate # cor T0, S1 | A
condat = pcor(cbind(ST[,2], A, ST[,3], A^2))$estimate # cor T0, T1 | A

condat[1,3] - condat0[1,3]/condat1[1,3] # CI conditional on A

condxt1 = pcor(cbind(ST[,3], X, ST[,1]))$estimate # cor T1, S1 | X
condxt0 = pcor(cbind(ST[,2], X, ST[,1]))$estimate # cor T0, S1 | X
condxt = pcor(cbind(ST[,2], X, ST[,3]))$estimate # cor T0, T1 | X

condxt[1,3] - condxt0[1,3]/condxt1[1,3] # CI conditional on X

cor(ST[,3], ST[,2]) - cor(ST[,1], ST[,2])/cor(ST[,1], ST[,3])

setind <- as.vector(mapply(rep, 1:iter, n)); set <- as.numeric(array_id);  X = X[setind == set]
samp=allsamp[setind==set,]; trt=alltrt[setind==set]
X=X[1:n]; A=A[1:n]
samp = samp[1:n, ]; trt = c(rep(0, n/2), rep(1, n/2))

sobs<-c(samp[(n/2):n,2]); tobs<-c(samp[1:(n/2-1),1],samp[(n/2):n,2])

ST<-samp; 
fsig<-function(n,s,j){  return(-(n)*log(s)+(-0.5*j)			)}
fsiginv<-function(n,s,j){  return((n)*log(s)+(-0.5*j)  )}
fdelt<-function(n, R, j){  return(-(n/2)*log(det(R))+(-0.5*j)				)}
fdeltBeta<-function(n, R, j){
  p = 5; q = 6 # want p>q
  a = -0.4 # lower bound of beta
  b = 1 # upper bound of beta
  return( -(n/2) * log(det(R)) + (-0.5 * j) + (q - 1) * log(b - R[2,3]) + # beta
            (p - 1) * log(R[2,3] - a) ) # alpha
}
tau1 = tau2 = tau3 = 0.25
S = diag(c(tau1, tau2, tau3))
s1=diag(S)[1]; s2=diag(S)[2]; s3=diag(S)[3]; SIG0Smat<-diag(.01,1)

holdmu<-matrix(rep(0,3*SIM),3,SIM); holdmu1<-matrix(rep(0,3*SIM),3,SIM)

holdS<-array(rep(0,3*3*SIM),dim=c(3,3,SIM)); holdT<-array(rep(0,3*3*SIM),dim=c(3,3,SIM))
holdR<-array(rep(0,3*3*SIM),dim=c(3,3,SIM)); holdR[,,1]=R; holdR[1,1,] = 1; holdS[,,1]=S; holdT[,,1]=ginv(S)
holdC<-array(rep(0,3*3*SIM),dim=c(3,3,SIM))
holdphi1<-(rep(0,1*SIM)); holdphi2<-(rep(0,1*SIM))
holdphi5<-(rep(0,1*SIM)); holdphi8<-(rep(0,1*SIM)); holdphi62<-(rep(0,1*SIM));
holdphi3<-(rep(0,1*SIM)); holdphi6<-(rep(0,1*SIM)); holdphi9<-(rep(0,1*SIM)); 
holdphi4<-(rep(0,1*SIM)); holdphi7<-(rep(0,1*SIM)); holdphi92<-(rep(0,1*SIM)); 
thetaS=theta0S<-matrix(c(1,0),nrow=2); thetaT=theta0T<-matrix(c(rep(0,4)),nrow=4)

S1 = seq(min(ST[, 1]), max(ST[,1]), 0.1); prob_t1t0<-matrix(0,(SIM), length(S1))

tauSS0=tauSS1=tauST0=tauST1=c(1); Xmat=cbind(1,X,A)
sim=2

surrogprob = function(x, S1, gamma0, gamma1, rho10, rhot, rho11, sigs1, sigt0, sigt1){
  pnorm(x, gamma0 + gamma1*S1,  sd = sqrt(sigt0^2 - sigt0^2*rho10^2 + sigt1^2 - 
                                            sigt1^2*rho11^2 - 2*rhot*sigt0*sigt1 + 2*rho10*rho11*sigt0*sigt1), lower.tail = FALSE)}

slope = slopemarginalizefull = intmarginalizefull = intmarginalizeclosed = slopemarginalizeclosed =
slopemarginalize = slopemarginalizewtcl = intmarginalizewtcl = 
intmarginalize = int = int1 = array(0,c((sim-1),1))

slope[1] = (holdR[1,3,1]*holdS[3,3,1]-holdR[1,2,1]*holdS[2,2,1])/(holdS[1,1,1])
slope[1] = (holdR[1,3,1]*(holdS[3,3,1] + tau)-holdR[1,2,1]*(holdS[2,2,1] + tau))/(holdS[1,1,1])
holdmu[,1] = c(phi1, phi4, phi7)
(holdmu[3,1]-holdmu[2,1])-slope[1]*(holdmu[1,1])
(holdmu1[3,1]-holdmu1[2,1])-slope[1]*(holdmu1[1,1])
slope[1]

slope[1] = (cor(ST[,3], ST[,1])*sd(ST[,3]) - cor(ST[,2], ST[,1])*sd(ST[,2]))/sd(ST[,1])
mean(ST[,3]) - mean(ST[,2]) - 
  (cor(ST[,3], ST[,1])*sd(ST[,3]) - cor(ST[,2], ST[,1])*sd(ST[,2]))/sd(ST[,1])*mean(ST[,1])
Age = A


# treatment effect at observed level
mean(ST[trt==1,3]) - mean(ST[trt==0,2])
mean(ST[trt==1,3] - X[trt==1]) - mean(ST[trt==0,2] - X[trt==0])
coef(lm(ST[,3]*trt + ST[,2]*(1-trt) ~ X + trt))[3]
coef(lm(ST[,3]*trt + ST[,2]*(1-trt) ~ X + trt + A + I(A^2)))[3]
coef(lm(ST[,3]*trt + ST[,2]*(1-trt) ~ trt + A + I(A^2)))[2]

summary(A)

holdphi1[1] = mean(ST[trt==1,1])
holdphi4[1] = as.numeric(coef(lm(ST[trt==0,2] ~ X[trt==0] + Age[trt==0] + I(Age[trt==0]^2)))[1])
holdphi5[1] = coef(lm(ST[trt==0,2] ~ X[trt==0] + Age[trt==0] + I(Age[trt==0]^2)))[2]
holdphi6[1] = coef(lm(ST[trt==0,2] ~ X[trt==0] + Age[trt==0] + I(Age[trt==0]^2)))[3]
holdphi7[1] = coef(lm(ST[trt==1,3] ~ X[trt==1] + Age[trt==1] + I(Age[trt==0]^2)))[1]
holdphi8[1] = coef(lm(ST[trt==1,3] ~ X[trt==1] + Age[trt==1] + I(Age[trt==0]^2)))[2]
holdphi9[1] = coef(lm(ST[trt==1,3] ~ X[trt==1] + Age[trt==1] + I(Age[trt==0]^2)))[3]

# surrogacy at smaller n level
lm(ST[,3] - ST[,2] ~ ST[,1])
lm(ST[,3] - ST[,2] ~ Age + ST[,1])
lm(ST[,3] - ST[,2] ~ Age + I(Age^2) + ST[,1])
lm(ST[,3] - ST[,2] ~ X + ST[,1])
lm(ST[,3] - ST[,2] ~ X + Age + ST[,1])
lm(ST[,3] - ST[,2] ~ X + Age + I(Age^2) + ST[,1])

# partial correlations
# is the correlation between the residuals
# are MVN, the partial correlation equals the conditional correlation.

condaxt1 = pcor(cbind(ST[,3], A, X, ST[,1], A^2))$estimate # cor T1, S1 | A, X
condaxt0 = pcor(cbind(ST[,2], A, X, ST[,1], A^2))$estimate # cor T0, S1 | A, X
condaxt = pcor(cbind(ST[,2], A, X, ST[,3], A^2))$estimate # cor T0, T1 | A, X

condaxt[1,4] - condaxt0[1,4]/condaxt1[1,4] # CI conditional on A, X

condaxt1 = pcor(cbind(ST[,3], A, X, ST[,1]))$estimate # cor T1, S1 | A, X
condaxt0 = pcor(cbind(ST[,2], A, X, ST[,1]))$estimate # cor T0, S1 | A, X
condaxt = pcor(cbind(ST[,2], A, X, ST[,3]))$estimate # cor T0, T1 | A, X

condaxt[1,4] - condaxt0[1,4]/condaxt1[1,4] # CI conditional on A, X

condat1 = pcor(cbind(ST[,3], A, ST[,1]))$estimate # cor T1, S1 | A
condat0 = pcor(cbind(ST[,2], A, ST[,1]))$estimate # cor T0, S1 | A
condat = pcor(cbind(ST[,2], A, ST[,3]))$estimate # cor T0, T1 | A

condat[1,3] - condat0[1,3]/condat1[1,3] # CI conditional on A

condxt1 = pcor(cbind(ST[,3], X, ST[,1]))$estimate # cor T1, S1 | X
condxt0 = pcor(cbind(ST[,2], X, ST[,1]))$estimate # cor T0, S1 | X
condxt = pcor(cbind(ST[,2], X, ST[,3]))$estimate # cor T0, T1 | X

condxt[1,3] - condxt0[1,3]/condxt1[1,3] # CI conditional on A

cor(ST[,3], ST[,2]) - cor(ST[,1], ST[,2])/cor(ST[,1], ST[,3])

lm(ST ~ X + Age + I(Age^2))

lm(ST[trt==1,] ~ X[trt==1] + Age[trt==1] + I(Age[trt==1]^2))
lm(ST[trt==0,] ~ X[trt==0] + Age[trt==0] + I(Age[trt==0]^2))

lm(ST[trt==1,1] ~ X[trt==1] + Age[trt==1] + I(Age[trt==1]^2))
lm(ST[trt==0,1] ~ X[trt==0] + Age[trt==0] + I(Age[trt==0]^2))

colMeans(ST[trt==1,])
colMeans(ST[trt==0,])
mean(Age[trt==1])
mean(Age[trt==0])
mean(X[trt==1])
mean(X[trt==0])

#ST[1:(n/2),1]<-0; ST[(n/2+1):n,2]<-0; ST[1:(n/2),3]<-0

while(sim<=SIM){
  
  mu = cbind(holdphi1[sim-1] + holdphi2[sim-1] * X + holdphi3[sim-1] * Age,
           as.numeric(holdphi4[sim-1]) + holdphi5[sim-1] * X + holdphi6[sim-1] * Age + holdphi62[sim-1] * Age^2,
           holdphi7[sim-1] + holdphi8[sim-1] * X + holdphi9[sim-1] * Age  + holdphi92[sim-1] * Age^2)
  
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
  var(ST[trt==1,])
  
  if(TRUE & sim == 2){
  ST = samp
  if(difference == TRUE) ST[,2:3] = ST[,2:3] - X
  }

  XmatS=c(rep(1, n))
  theta0S = 0
  
  #### Estimate S1, T0, T1|X
  Xmat = rep(1, n) # n by k
  Sigmats = Sig # 3 x 3
  
  if(Xfit == FALSE){
    
    Lambda0t = matrix(c(rep(0, 1)), nrow = 1); diag(Lambda0t) = .1
    Lambdan <- t(Xmat)%*%Xmat + Lambda0t # put zeroes in
    Lambdaninv = round(ginv(Lambdan), 9)
    B0 = matrix(c(rep(0, 3)), nrow = 1)
    Lambdans = t(Xmat)%*%Xmat + Lambda0t
    betaTn = ginv(t(Xmat)%*%Xmat + Lambda0t) %*% (t(Xmat)%*% ST[,1:3] + Lambda0t %*% B0)
    betaT = rmatrixnorm(betaTn, Lambdaninv, as.symmetric.matrix(Sig))      # dim nxk (3 x 2), nxn (3x3), kxk (2x2)         
    
    holdphi1[sim]=betaT[1] 
    
    holdphi4[sim]=betaT[2]
    
    holdphi7[sim]=betaT[3]
    
    XmatS = Xmat; SIG0Smat = c(0.1)
    v<-ginv(SIG0Smat+as.numeric(tauST0)*(t(XmatS)%*%XmatS))
    m<-v%*%(tauST0*t(XmatS)%*%ST[,2])  # more generally, diag(t0mat)%*%m0
    betaT<-c(rmvnorm(1,m,v))
    holdphi4[sim]=betaT[1]
    tmp2<-(ST[,2])-XmatS*holdphi4[sim]
    
    v<-ginv(SIG0Smat+as.numeric(tauST1)*(t(XmatS)%*%XmatS))
    m<-v%*%(tauST1*t(XmatS)%*%ST[,3])  # more generally, diag(t0mat)%*%m0
    betaT<-c(rmvnorm(1,m,v))
    holdphi7[sim]=betaT[1]
    tmp3<-(ST[,3])-XmatS*holdphi7[sim]
    
    v<-ginv(SIG0Smat+as.numeric(tauSS0)*(t(XmatS)%*%XmatS))
    m<-v%*%(tauSS0*t(XmatS)%*%ST[,1])  # more generally, diag(t0mat)%*%m0
    betaT<-c(rmvnorm(1,m,v))
    holdphi1[sim]=betaT[1]
    tmp1<-(ST[,1])-XmatS*holdphi1[sim]
    
    
    
  }
  
  if(Xfit == TRUE){
    Xmat = cbind(1, X, Age, Age^2) # n by k
    
    Lambda0t = matrix(c(rep(0, 16)), nrow = 4); diag(Lambda0t) = .1
    Lambdan <- t(Xmat)%*%Xmat + Lambda0t # put zeroes in
    Lambdaninv = round(ginv(Lambdan), 9)
    Lambdans = t(Xmat)%*%Xmat + Lambda0t
    betaTn = ginv(t(Xmat)%*%Xmat + Lambda0t) %*% (t(Xmat)%*% ST)
    betaT = rmatrixnorm(betaTn, Lambdaninv, as.symmetric.matrix(Sig))      # dim nxk (3 x 2), nxn (3x3), kxk (2x2)         
    
    betaTn = ginv(t(Xmat)%*%Xmat + Lambda0t) %*% (t(Xmat)%*% ST[,2:3])
    betaT = rmatrixnorm(betaTn, Lambdaninv, as.symmetric.matrix(Sig[2:3, 2:3]))      # dim nxk (3 x 2), nxn (3x3), kxk (2x2)         
  
    #holdphi1[sim]=betaT[1,1]
    #holdphi2[sim]=betaT[2,1] = 0
    #holdphi3[sim]=betaT[3,1] = 0
    
    holdphi4[sim]=betaT[1,2]
    holdphi5[sim]=betaT[2,2]
    holdphi6[sim]=betaT[3,2]
    holdphi62[sim]=betaT[4,2]
    
    #holdphi7[sim]=betaT[1,3]
    #holdphi8[sim]=betaT[2,3]
    #holdphi9[sim]=betaT[3,3]
    #holdphi92[sim]=betaT[4,3]
    
    holdphi4[sim]=betaT[1,1]
    holdphi5[sim]=betaT[2,1]
    holdphi6[sim]=betaT[3,1]
    holdphi62[sim]=betaT[4,1]
    
    holdphi7[sim]=betaT[1,2]
    holdphi8[sim]=betaT[2,2]
    holdphi9[sim]=betaT[3,2]
    holdphi92[sim]=betaT[4,2]
    
    v<-ginv(SIG0Smat+as.numeric(tauSS0)*(t(XmatS)%*%XmatS))
    m<-v%*%(tauSS0*t(XmatS)%*%ST[,1])  # more generally, diag(t0mat)%*%m0
    holdphi1[sim] = c(rmvnorm(1,m,v))
  }
  
  tmp1<-(ST[,1])-XmatS*holdphi1[sim]
  tmp2<-(ST[,2])-XmatS*holdphi4[sim]
  tmp3<-(ST[,3])-XmatS*holdphi7[sim]  
  
  if(Xfit == TRUE){
    #tmp1<-(ST[,1])-Xmat%*%betaT[,1]
    tmp2<-(ST[,2])-Xmat%*%betaT[,1]
    tmp3<-(ST[,3])-Xmat%*%betaT[,2]
  }
  resid<-cbind((tmp1),(tmp2),(tmp3))
  
  ####calc s's
  s1 = rinvgamma(1, shape = a + n/2, scale = (sum(tmp1^2)/2 + b))  
  s2 = rinvgamma(1, shape = a + n/2, scale = (sum(tmp2^2)/2 + b))  
  s3 = rinvgamma(1, shape = a + n/2, scale = (sum(tmp3^2)/2 + b))  
  
  S = diag(c(sqrt(s1), sqrt(s2), sqrt(s3)))
  holdS[,,sim] = S
  
  ###r13 rho11
  R = holdR[,,sim-1]
  diag(R) = 1
  r23 = holdR[2,3,sim-1]; r13 = holdR[1,3,sim]; r12 = holdR[1,2,sim-1]
  
  U13 = (-2*r12*r23 - sqrt(4*r12^2*r23^2 + 4*(1 - r23^2 - r12^2)))/(-2)
  L13 = (-2*r12*r23 + sqrt(4*r12^2*r23^2 + 4*(1 - r23^2 - r12^2)))/(-2)
  if(TRUE){
    if(is.na(L13)|is.na(U13)) next
    
    low13=ceiling(50*max(min(L13,U13))); up13=min(floor(50*max(L13,U13)))
    d13=up13-low13+1
    
    fr13<-matrix(rep(0,d13*4),d13,4)
    for (k in low13:up13){
      r<-  k/50; R[1,3]=r; R[3,1]=r; fr13[(k-low13+1),1]=r
      # R[1,2] = R[2,1] = r
      if(condindfit == TRUE) R[1,2] = R[2,1] = r*R[2,3]
      
      if(any(is.nan(eigen(R)$values))){ next}
      if(any(eigen(R)$values<0)) next
      summand = apply(resid, 1, function(resid) resid%*%ginv(S%*%R%*%S)%*%resid)
      
      j3 = sum(summand)
      
      fr13[(k-low13+1),2]<-fdelt(n, R, j3)
    }
    
    fr13<-fr13[!is.na(fr13[,2]),]
    fr13<-fr13[!is.infinite(fr13[,2]),]
    fr13<-fr13[complete.cases(fr13),]
    fr13[,2]=fr13[,2]-median(fr13[,2])
    fr13[,2]=exp(fr13[,2])
    fr13<-fr13[!is.infinite(fr13[,2]),]
    m<-which(fr13[,2]==max(fr13[,2]))[1]
    fr13[m,]
    x1=max(1, (m-10)); x2=min(length(fr13[,1]), (m+10))
    
    low13=500*fr13[x1,1]; up13=500*fr13[x2,1]
    d13=up13-low13+1
    
    fr13<-matrix(rep(0,d13*4),d13,4)
    for (k in low13:up13){
      r<-  k/500; R[1,3]=r; R[3,1]=r; fr13[(k-low13+1),1]=r
      if(condindfit == TRUE) R[1,2] = R[2,1] = r*R[2,3]
      if(any(is.nan(eigen(R)$values))){ next}
      if(any(eigen(R)$values<0)) next
      summand = apply(resid, 1, function(resid) resid%*%ginv(S%*%R%*%S)%*%resid)
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
  }
  R[1,3]=r13; R[3,1]=r13
  R[1,2] = R[2,1] = r13*R[2,3]
  
  ### r23 rt
  if(TRUE){
    r23 = R[2,3]; r13 = R[1,3]; r12 = R[1,2]
    U23 = (-2*r12*r13-sqrt(4*r12^2*r13^2 + 4*(1 - r13^2 - r12^2)))/(-2)
    L23 = (-2*r12*r13+sqrt(4*r12^2*r13^2 + 4*(1 - r13^2 - r12^2)))/(-2)
    if(is.na(L23)|is.na(U23)) next
    
    low23=ceiling(50*max(min(L23,U23), -0.4)); up23=min(floor(50*max(L23,U23)))
    d23=up23-low23+1; fr23<-matrix(rep(0,d23*4),d23,4)
    for (k in low23:up23){
      r<-  k/50; R[2,3]=r; R[3,2]=r; fr23[(k-low23+1),1]=r
      
      if(condindfit == TRUE) R[1,2] = R[2,1] = r*R[1,3]
      if(any(is.nan(eigen(R)$values))){ next}
      if(any(eigen(R)$values<0)) next
      summand = apply(resid, 1, function(resid) resid%*%ginv(S%*%R%*%S)%*%resid)
      j3 = sum(summand)
      fr23[(k-low23+1),2]<-fdelt(n, R, j3)
    }
    
    fr23<-fr23[!is.na(fr23[,2]),]
    fr23<-fr23[!is.infinite(fr23[,2]),] 
    fr23<-fr23[complete.cases(fr23),]
    fr23[,2]=fr23[,2]-median(fr23[,2])
    fr23[,2]=exp(fr23[,2])
    fr23<-fr23[!is.infinite(fr23[,2]),]  
    m<-which(fr23[,2]==max(fr23[,2]))[1]
    fr23[m,]; plot(fr23[,1], fr23[,2])
    x1=max(1, (m-10)); x2=min(length(fr23[,1]), (m+10))
    low23=500*fr23[x1,1]; up23=500*fr23[x2,1]
    d23=up23-low23+1
    
    fr23<-matrix(rep(0,d23*4),d23,4)
    for (k in low23:up23){
      r <-  k/500; R[2,3]=r; R[3,2]=r; fr23[(k-low23+1),1] = r
      if(condindfit == TRUE) R[1,2] = R[2,1] = r*R[1,3]
      if(any(is.nan(eigen(R)$values))){ next}
      if(any(eigen(R)$values<0)) next      
      summand = apply(resid, 1, function(resid) resid%*%ginv(S%*%R%*%S)%*%resid)
      
      j3 = sum(summand)
      
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
  }
  R[2,3]=r23; R[3,2]=r23
  
  R[1,2] = R[2,1] = R[1,3] * R[2,3]
  
  if(condindfit == FALSE){
    U12 = (-2*r13*r23-sqrt(4*r13^2*r23^2 + 4*(1 - r23^2 - r13^2)))/(-2)
    L12 = (-2*r13*r23+sqrt(4*r13^2*r23^2 + 4*(1 - r23^2 - r13^2)))/(-2)

    det(R[1:2, 1:2])-det(R[2:3, 2:3])+det(R[c(1,3), c(1,3)])
    det(R[1:2, 1:2])+det(R[2:3, 2:3])+det(R[c(1,3), c(1,3)])
    
    if(is.na(L12)|is.na(U12)) next
    
    low12=ceiling(50*max(min(L12,U12), -1)); #up12=(floor(100*max(L12,min(U12))))
    up12=floor(min(50*max(L12,min(U12),1)))
    d12=up12-low12+1
    
    fr12<-matrix(rep(0,d12*4),d12,4)
    for (k in low12:up12){
      r<- k/50;  R[1,2]=r; R[2,1]=r
      fr12[(k-low12+1),1]=r
      summand = apply(resid, 1, function(resid) resid%*%ginv(S%*%R%*%S)%*%resid)
      
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
    low12=500*fr12[x1,1]; up12=500*fr12[x2,1]
    d12=up12-low12+1
    
    fr12<-matrix(rep(0,d12*4),d12,4)
    for (k in low12:up12){
      r<- k/500;  R[1,2]=r;  R[2,1]=r
      fr12[(k-low12+1),1]=r
      summand = apply(resid, 1, function(resid) resid%*%ginv(S%*%R%*%S)%*%resid)
      
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
    
    tauST1<-holdS[3,3,sim-1]^2
    tauST0<-holdS[2,2,sim-1]^2
    tauSS0<-holdS[1,1,sim-1]^2
    Xmat = cbind(1, X, Age, Age^2) # n by k
    
    v<-ginv(SIG0Smat[1,1]+as.numeric(ginv(tauST0))*(t(Xmat[trt==0,])%*%Xmat[trt==0,]))
    
    m<-v%*%(as.numeric(ginv(tauST0))*t(Xmat[trt==0,])%*%ST[trt == 0,2])  # more generally, diag(t0mat)%*%m0
    betaT<-c(rmvnorm(1,m,v^2))
    
    #betaT = as.numeric(coef(lm(ST[trt==0, 2] ~ X[trt==0] + Age[trt==0] + I(Age[trt==0]^2))))

    holdphi4[sim]=betaT[1]
    holdphi5[sim]=betaT[2]
    holdphi6[sim]=betaT[3]
    holdphi62[sim]=betaT[4]
    
    tmp2<-(ST[trt == 0,2])-Xmat[trt==0,]%*%betaT
    
    v<-ginv(SIG0Smat[1,1]+as.numeric(ginv(tauST1))*(t(Xmat[trt==1,])%*%Xmat[trt==1,]))

    m<-v%*%(as.numeric(ginv(tauST1))*t(Xmat[trt==1,])%*%ST[trt == 1,3])  # more generally, diag(t0mat)%*%m0
    betaT<-c(rmvnorm(1,m,v^2))
    
    #betaT = as.numeric(coef(lm(ST[trt==1, 3] ~ X[trt==1] + Age[trt==1] + I(Age[trt==1]^2))))
    
    holdphi7[sim]=betaT[1]
    holdphi7[sim] = as.numeric(coef(lm(ST[trt==1, 3] ~ X[trt==1] + Age[trt==1] + I(Age[trt==1]^2))))[1]
    holdphi8[sim]=betaT[2]
    holdphi9[sim]=betaT[3]
    holdphi92[sim]=betaT[4]
    tmp3<-(ST[trt == 1,3])-Xmat[trt==1,]%*%betaT
    
    v<-ginv(SIG0Smat[1,1]+as.numeric(tauSS0)*(t(Xmat[trt==1,1])%*%Xmat[trt==1,1]))
    m<-v%*%(tauSS0*t(Xmat[trt==1,1])%*%ST[trt == 1,1])  # more generally, diag(t0mat)%*%m0
    betaT<-c(rmvnorm(1,m,v^2))
    holdphi1[sim]=betaT[1]
    tmp1<-(ST[trt == 1,1])-Xmat[trt==1,1] * holdphi1[sim]
    
    if(FALSE){ 
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
        
        
      }
    }
    
    holdmu[,sim] = c(holdphi1[sim], holdphi4[sim], holdphi7[sim])
    
    # inv gamma priors for S
    s1 =   rinvgamma(1, shape = a + n/4, scale = (sum(tmp1^2)/2 + b))
    s2 =   rinvgamma(1, shape = a + n/4, scale = (sum(tmp2^2)/2 + b))  
    s3 =   rinvgamma(1, shape = a + n/4, scale = (sum(tmp3^2)/2 + b))  
    
    holdS[,,sim] = S = diag(c(sqrt(s1), sqrt(s2), sqrt(s3)))
    
    ###r13 rho11
    R = holdR[,,sim-1]
    diag(R) = 1
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
    
    #R[2,3] = R[3,2] = r23 = holdR[2,3,sim] = rBeta_ab(1, shape1 = 5, shape2 = 6, a = -0.4, b = 1)
    #R[2,3] = R[3,2] = r23 = holdR[2,3,sim] = 0.22
    #R[2,3] = R[3,2] = r23 = holdR[2,3,sim] = 0.04
    
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
  
  if(any(is.nan(eigen(R)$values))){ next}
  if(any(eigen(R)$values<0)) next
  if(any(R>1.000000001)) next
  Sig = S%*%R%*%S

  holdmu[,sim] = c(holdphi1[sim], holdphi4[sim], holdphi7[sim])
  holdR[,,sim] = R
  
  slope[sim]=(holdR[1,3,sim]*holdS[3,3,sim]-holdR[1,2,sim]*holdS[2,2,sim])/(holdS[1,1,sim])
  slope[1]=(holdR[1,3,1]*holdS[3,3,1]-holdR[1,2,1]*holdS[2,2,1])/(holdS[1,1,1])
  
  if(Xfit == TRUE){ # decide how we want to do this part -> age at 4 years
    holdmu1[,sim] = c(holdphi1[sim] + holdphi2[sim], holdphi4[sim] + holdphi5[sim], holdphi7[sim] + holdphi8[sim])
    holdmu1[,sim] = c(holdphi1[sim], holdphi4[sim] + 5*holdphi6[sim] + 5^2*holdphi62[sim], holdphi7[sim] + 5*holdphi9[sim] +
                        5^2*holdphi92[sim])
    
    holdmu1[,sim] = c(holdphi1[sim], holdphi4[sim] + 5*holdphi6[sim] + 5^2*holdphi62[sim] + 22*holdphi5[sim], 
                      holdphi7[sim] + 5*holdphi9[sim] + 5^2*holdphi92[sim] + 22*holdphi8[sim])
    
    int1[sim]=(holdmu1[3,sim]-holdmu1[2,sim])-slope[sim]*(holdmu1[1,sim])

    eta1h = mean(ST[,1]); eta2h = mean(ST[,2]); eta3h = mean(ST[,3]); eta4h = mean(X); eta5h = mean(A) # mean terms
    tau1h = sd(ST[,1]); tau2h = sd(ST[,2]); tau3h = sd(ST[,3]); tau4h = sd(X); tau5h = sd(A) # diagonal sd terms
    pi10h = cor(ST[,1], ST[,2]); pi11h = cor(ST[,1], ST[,3]); pi1Xh = cor(ST[,1], X)
    piTh = cor(ST[,2], ST[,3]); piX0h = cor(ST[,2], X)
    piX1h = cor(ST[,3], X); piA1h = cor(ST[,3], Age); piA0h = cor(ST[,2], Age); pi1Ah = cor(ST[,1], Age)
    
    ### make conditional on age
    eta2h = coef(lm(ST[,2] ~ Age + I(Age^2)))[1]; eta3h = coef(lm(ST[,2] ~ Age + I(Age^2)))[1]; eta4h = mean(X); eta5h = mean(A) # mean terms
    tau1h = sd(ST[,1]); tau2h = S[2,2]; tau3h = S[3,3]; tau4h = sd(X); tau5h = sd(A) # diagonal sd terms
    pi10h = pcor(cbind(ST[,2], Age, ST[,1], Age^2))$estimate[1,3]; pi11h = pcor(cbind(ST[,3], Age, ST[,1], Age^2))$estimate[1,3]; 
    pi1Xh = cor(ST[,1], X)
    piTh = pcor(cbind(ST[,2], Age, ST[,3], Age^2))$estimate[1,3]; piX0h = pcor(cbind(ST[,2], Age, X, Age^2))$estimate[1,3]
    piX1h = pcor(cbind(ST[,3], Age, X, Age^2))$estimate[1,3]; 
    
    s =  seq(min(samp[,1], -1), max(samp[,1]), length.out = n) # grid should stay same every time
    intgrid = probgrid = probgrid2 = probgridint = weightgrid = weightgridclosed = Xgrid = sgrid = matrix(s, length(s), length(s))
    Xgrid[,1] = X 
    for(l in 2:ncol(sgrid)){    Xgrid[,l] = c(X[l:(n)], X[1:(l - 1)])   }
    
    for(l in 1:ncol(sgrid)){
      
      Xpermute = Xgrid[,l]
      # needs to be |Age, Age^2
      
      intgrid[,l] = intx = ((coef(lm(ST[,3] ~ X + Age + I(Age^2)))[1] + coef(lm(ST[,3] ~ X + Age + I(Age^2)))[2]*Xpermute) - 
                              (coef(lm(ST[,2] ~ X + Age + I(Age^2)))[1] + coef(lm(ST[,2] ~ X + Age + I(Age^2)))[2]*Xpermute))-slope[sim] * (
                                holdphi1[sim])
      sumSX = 0*Xpermute; xindex = 1
      
      for(xindex in 1:nrow(sgrid)){                       
        sumSX[xindex] = exp(-1/(2 * S[1,1]^2) * (sgrid[xindex, l] - (c(1, Xgrid[xindex, l]) %*% c(holdphi1[sim], holdphi2[sim]) ) )^2) # gives correct answer
      } # gives same as closed form

      weightgrid[,l] = sumSX  # /sum(sumSX)
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
    if(sim <= 3){
      plot(s, rowSums(probgrid)); lines(s, int[sim] + slope[sim]*s, col = "red"); if(Xfit == TRUE){lines(s, int1[sim] + slope[sim]*s, col = "blue")}
      plot(s, rowSums(probgrid2)); lines(s, int[sim] + slope[sim]*s, col = "red"); if(Xfit == TRUE){lines(s, int1[sim] + slope[sim]*s, col = "blue")}
    }
    
    intmarginalize[sim] = lm((rowSums(probgrid)) ~ s)$coef[1]
    
    slopemarginalize[sim] = lm((rowSums(probgrid)) ~ s)$coef[2]
    
    intmarginalizewtcl[sim] = lm((rowSums(probgrid2)) ~ s)$coef[1]
    
    slopemarginalizewtcl[sim] = lm((rowSums(probgrid2)) ~ s)$coef[2]
    
    # closed form version
    A = holdphi7[sim] - holdphi4[sim] - slope[sim]*holdphi1[sim]
    B = holdphi8[sim] - holdphi5[sim] - slope[sim]*holdphi2[sim]
    
    gamma0star = A + B*(eta4h - tau1h^(-1)*tau4h*pi1Xh*eta1h) # aligns with previous belief
    
    gamma1star = B*(tau1h^(-1)*tau4h*pi1Xh) + slope[sim]
    if(sim <= 3){
      plot(s,  gamma0star + gamma1star*s)
      lines(s, int[sim] + slope[sim]*s, col = "red");
      if(Xfit == TRUE) lines(s, int1[sim] + slope[sim]*s, col = "blue");
    }
    
    intmarginalizeclosed[sim] = gamma0star
    slopemarginalizeclosed[sim] = gamma1star
    
    ### repeat to integrate out A
    delta1h = mean(ST[,1])
    rhos1xh = cor(ST[,1], X)
    rhos1ah = cor(ST[,1], Age)
    tau5h = sd(Age)
    B = holdphi7[sim] - holdphi4[sim] - slope[sim] * holdphi1[sim]
    M = holdphi8[sim] - holdphi5[sim] - slope[sim] * holdphi2[sim]
    N = holdphi9[sim] - holdphi6[sim] - slope[sim] * holdphi3[sim]
    intmarginalizefull[sim] = B + M*(eta4h + holdS[1,1,sim]^(-1)*tau4h*rhos1xh*delta1h) + N*
      (eta5h + holdS[1,1,sim]^(-1)*tau5h*rhos1ah*delta1h)
    
    slopemarginalizefull[sim] = -M*holdS[1,1,sim]^(-1)*tau4h*rhos1xh - N*
      holdS[1,1,sim]^(-1)*tau5h*rhos1ah + slope[sim]
    
  }
  
  int[sim]=(holdmu[3,sim]-holdmu[2,sim])-slope[sim]*(holdmu[1,sim])
  int1[sim]
  
  prob_t1t0[sim,] = surrogprob(0, S1, gamma0 = int[sim], gamma1 = slope[sim], 
                               rho10 = R[1,2], rhot = R[2,3], rho11 = R[1,3], sigs1 = holdS[1,1,sim],
                               sigt0 = holdS[2,2,sim], sigt1 = holdS[3,3,sim])
  
  
  slope[1] = (holdR[1,3,1]*holdS[3,3,1]-holdR[1,2,1]*holdS[2,2,1])/(holdS[1,1,1])
  R[1,2] = R[2,1] = R[1,3] * R[2,3]
  
  if(any(eigen(R)$values<0)) next; if(any(abs(R)>1)) next
  holdR[,, sim] = R
  
  print(R)
  print(S)
  print(slope[sim])
  print(slopemarginalizeclosed[sim])
  print(slopemarginalizefull[sim])
  
  print(int[sim])
  print(int1[sim])
  
  print(intmarginalizeclosed[sim])
  print(intmarginalizefull[sim])
  
  sim=sim+1
  print(sim)
}    

mean(slopemarginalizefull)
mean(intmarginalizefull)
mean(slopemarginalizeclosed)
mean(intmarginalizeclosed)
mean(slopemarginalize)
mean(intmarginalize)

mean(slope)
mean(int)
mean(int1)

slope[1]=(holdR[1,3,1]*holdS[3,3,1]-holdR[1,2,1]*holdS[2,2,1])/(holdS[1,1,1])
int[1]=(holdmu[3,1]-holdmu[2,1])-slope[1]*(holdmu[1,1])


prob_t1t0[1,] = surrogprob(0, S1, gamma0 = int[1], gamma1 = slope[1], 
                           rho10 = holdR[1,2,1], rhot = holdR[2,3,1], rho11 = holdR[1,3,1], sigs1 = holdS[1,1,sim],
                           sigt0 = holdS[2,2,1], sigt1 = holdS[3,3,1])

S1 = seq(-max(ST[, 1]), max(ST[,1]), 0.1)

prob_t1t0l = apply(X = prob_t1t0, MARGIN = 2, FUN = quantile, probs = 0.025, na.rm=T)
prob_t1t0u = apply(X = prob_t1t0, MARGIN = 2, FUN = quantile, probs = 0.975, na.rm=T)

#plot(S1, colMeans(prob_t1t0, na.rm = T), xlab = "S1 Value", ylab = "Probability(T1-T0>0|S1)", main = 
#       "Treatment Efficacy by Surrogate", ylim = c(0, 1)); lines(S1, prob_t1t0l); lines(S1, prob_t1t0u)


intl = quantile(int[burnin:sim-1],  probs = 0.025, na.rm=T) + quantile(slope[burnin:sim-1],  probs = 0.025, na.rm=T)*S1
intu = quantile(int[burnin:sim-1],  probs = 0.975, na.rm=T) + quantile(slope[burnin:sim-1],  probs = 0.975, na.rm=T)*S1

intl = sapply(X = S1,  FUN = quantile, probs = 0.025, na.rm=T, int + slope*X)
intu = sapply(X = S1,  FUN = quantile, probs = 0.975, na.rm=T, int + slope*X)

ints = matrix(0,nrow = (length(int)), ncol = length(S1))
for(i in 1:length(int)){
  ints[i,] = int[i] + slope[i]*S1
}

intl = quantile(ints,  probs = 0.025, na.rm=T)*S1
intu = quantile(ints,  probs = 0.975, na.rm=T)*S1

S1 = seq(0, max(ST[,1]), 0.1)

prob_t1t0 = surrogprob(0, S1, gamma0 = int[1], gamma1 = slope[1], 
                       rho10 = holdR[1,2,1], rhot = holdR[2,3,1], rho11 = holdR[1,3,1], sigs1 = holdS[1,1,sim],
                       sigt0 = holdS[2,2,1], sigt1 = holdS[3,3,1])

plot(S1, prob_t1t0, xlab = "S1 Value", ylab = "Probability(T1-T0>0|S1)", main = 
       "Treatment Efficacy by Surrogate,\nTruth", ylim = c(0, 1))


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
                   p01l=numeric(1),p01u=numeric(1),p10l=numeric(1),p10u=numeric(1),
                   r1=numeric(1), SEr1=numeric(1), r2=numeric(1), SEr2=numeric(1),r3=numeric(1),SEr3=numeric(1),
                   r1l=numeric(1), r1u=numeric(1),r2l=numeric(1), r2u=numeric(1),
                   r3l=numeric(1),r3u=numeric(1), b0=numeric(1), b0SE=numeric(1), b1=numeric(1), b1SE=numeric(1))

params[1]<-sqrt(var(holdphi4[burnin:(sim-1)]))
params[2]<-sqrt(var(holdphi7[burnin:(sim-1)]))
params[3]<-sqrt(var(holdphi9[burnin:(sim-1)]))
params[4]<-sqrt(var(holdphi1[burnin:(sim-1)]))
params[5]<-sqrt(var(holdphi2[burnin:(sim-1)]))
params[6]<-sqrt(var(holdphi5[burnin:(sim-1)]))
params[7]<-sqrt(var(holdphi6[burnin:(sim-1)]))
params[8]<-sqrt(var(holdphi8[burnin:(sim-1)]))

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

#params[45]<-mean(holdR2[1,2,burnin:(sim-1)])
#params[46]<-sqrt(var(holdR2[1,2,burnin:(sim-1)]))
#params[47]<-mean(holdR2[1,3,burnin:(sim-1)])
#params[48]<-sqrt(var(holdR2[1,3,burnin:(sim-1)]))
#params[49]<-mean(holdR2[2,3,burnin:(sim-1)])
#params[50]<-sqrt(var(holdR2[2,3,burnin:(sim-1)]))
#params[51]<-quantile(holdR2[1,2,burnin:(sim-1)],  probs = 0.025,na.rm=TRUE)
#params[52]<-quantile(holdR2[1,2,burnin:(sim-1)],  probs = 0.975,na.rm=TRUE)
#params[53]<-quantile(holdR2[1,3,burnin:(sim-1)],  probs = 0.025,na.rm=TRUE)
#params[54]<-quantile(holdR2[1,3,burnin:(sim-1)],  probs = 0.975,na.rm=TRUE)
#params[55]<-quantile(holdR2[2,3,burnin:(sim-1)],  probs = 0.025,na.rm=TRUE)
#params[56]<-quantile(holdR2[2,3,burnin:(sim-1)],  probs = 0.975,na.rm=TRUE)
#params[57]<-mean(holdB0[burnin:(sim-1)])
#params[58]<-sqrt(var(holdB0[burnin:(sim-1)]))
#params[59]<-mean(holdB1[burnin:(sim-1)])
#params[60]<-sqrt(var(holdB1[burnin:(sim-1)]))


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
               intmargfulllcl=numeric(1), intmargfullucl=numeric(1),
               slopemargfull_meancd=numeric(1), slopemargfull_SEcd=numeric(1), 
               slopemargfulllcd=numeric(1), slopemargfullucd=numeric(1),
               slopemargfull_mean=numeric(1), slopemargfull_SE=numeric(1))

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
PS[47]<-quantile(intmarginalizefull[burnin:sim-1],  probs = 0.025, na.rm=T)
PS[48]<-quantile(intmarginalizefull[burnin:sim-1],  probs = 0.975, na.rm=T)
PS[49]<-mean(intmarginalizefull[burnin:sim-1], na.rm=T)
PS[50]<-sqrt(var(intmarginalizefull[burnin:sim-1], na.rm=T))
PS[51]<-quantile(slopemarginalizefull[burnin:sim-1],  probs = 0.025, na.rm=T)
PS[52]<-quantile(slopemarginalizefull[burnin:sim-1],  probs = 0.975, na.rm=T)
PS[53]<-mean(slopemarginalizefull[burnin:sim-1], na.rm=T)
PS[54]<-sqrt(var(slopemarginalizefull[burnin:sim-1], na.rm=T))

estcoef<-data.frame(beta0mean=numeric(1), beta1mean=numeric(1),omega1mean=numeric(1), alpha0mean=numeric(1),
                    alpha1mean=numeric(1),psi1mean=numeric(1),psi2mean=numeric(1),omega2mean=numeric(1),
                    beta0l=numeric(1),beta0u=numeric(1), beta1l=numeric(1),beta1u=numeric(1),omega1l=numeric(1),
                    omega1u=numeric(1), alpha0l=numeric(1),alpha0u=numeric(1),
                    alpha1l=numeric(1),alpha1u=numeric(1),psi1l=numeric(1),psi1u=numeric(1),
                    psi2l=numeric(1),psi2u=numeric(1),omega2l=numeric(1),omega2u=numeric(1), 
                    beta0ind=numeric(1), beta1ind=numeric(1),omega1ind=numeric(1), alpha0ind=numeric(1),
                    alpha1ind=numeric(1),psi1ind=numeric(1),psi2ind=numeric(1),omega2ind=numeric(1))

estcoef[1]=mean(holdphi4[burnin:sim-1], na.rm=T)
estcoef[2]=mean(holdphi7[burnin:sim-1], na.rm=T)
estcoef[3]=mean(holdphi5[burnin:sim-1], na.rm=T)
estcoef[4]=mean(holdphi9[burnin:sim-1], na.rm=T)
estcoef[5]=mean(holdphi1[burnin:sim-1], na.rm=T)
estcoef[6]=mean(holdphi8[burnin:sim-1], na.rm=T)
estcoef[7]=mean(holdphi2[burnin:sim-1], na.rm=T)
estcoef[8]=omega2mean=mean(holdphi8[burnin:sim-1], na.rm=T)
estcoef[9]<-quantile(holdphi4[burnin:sim-1],  probs = 0.025, na.rm=T)
estcoef[10]<-quantile(holdphi4[burnin:sim-1],  probs = 0.975, na.rm=T)
estcoef[11]<-quantile(holdphi7[burnin:sim-1],  probs = 0.025, na.rm=T)
estcoef[12]<-quantile(holdphi7[burnin:sim-1],  probs = 0.975, na.rm=T)
estcoef[13]<-quantile(holdphi5[burnin:sim-1],  probs = 0.025, na.rm=T)
estcoef[14]<-quantile(holdphi5[burnin:sim-1],  probs = 0.975, na.rm=T)
estcoef[15]<-quantile(holdphi6[burnin:sim-1],  probs = 0.025, na.rm=T)
estcoef[16]<-quantile(holdphi6[burnin:sim-1],  probs = 0.975, na.rm=T)
estcoef[17]<-quantile(holdphi1[burnin:sim-1],  probs = 0.025, na.rm=T)
estcoef[18]<-quantile(holdphi1[burnin:sim-1],  probs = 0.975, na.rm=T)
estcoef[19]<-quantile(holdphi9[burnin:sim-1],  probs = 0.025, na.rm=T)
estcoef[20]<-quantile(holdphi9[burnin:sim-1],  probs = 0.975, na.rm=T)
estcoef[21]<-quantile(holdphi2[burnin:sim-1],  probs = 0.025, na.rm=T)
estcoef[22]<-quantile(holdphi2[burnin:sim-1],  probs = 0.975, na.rm=T)
estcoef[23]<-quantile(holdphi8[burnin:sim-1],  probs = 0.025, na.rm=T)
estcoef[24]<-quantile(holdphi8[burnin:sim-1],  probs = 0.975, na.rm=T)
estcoef[25]=as.numeric(estcoef[9]<holdmu[2,1]&estcoef[10]>holdmu[2,1])
estcoef[26]=as.numeric(estcoef[11]<holdmu[3,1]&estcoef[12]>holdmu[3,1])
estcoef[27]=as.numeric(estcoef[13]< holdphi5[1] &estcoef[14]> holdphi5[1])
estcoef[28]=as.numeric(estcoef[15]<holdmu[1,1]&estcoef[16]>holdmu[1,1])
estcoef[29]=as.numeric(estcoef[17]<holdmu[1,1]&estcoef[18]>holdmu[1,1])
#estcoef[30]=as.numeric(estcoef[19]<1&estcoef[20]>1)
estcoef[31]=as.numeric(estcoef[21]< holdphi2[1] &estcoef[22]> holdphi2[1])
estcoef[32]=as.numeric(estcoef[23]< holdphi8[1]&estcoef[24]>holdphi8[1])

fname <- paste('params','coefsimVaccine.n',n,array_id,'.txt',sep="")
write.table(params, file=fname, sep="\t", row.names=F, col.names=T)
fname2 <- paste('PS','coefsimVaccine.n',n,array_id,'.txt',sep="")
write.table(PS, file=fname2, sep="\t", row.names=F, col.names=T)
fname3 <- paste('prentice','coefsim4Vaccine.n',n,array_id,'.txt',sep="")
###write.table(prentice, file=fname3, sep="\t", row.names=F, col.names=T)
fname4 <- paste('postpren','coefsimVaccine.n',n, array_id,'.txt',sep="")
###write.table(pren.post, file=fname4, sep="\t", row.names=F, col.names=T)
fname5 <- paste('naivemodels','coefsimVaccine.n',n,array_id,'.txt',sep="")
###write.table(naiveresults, file=fname5, sep="\t", row.names=F, col.names=T)
fname6 <- paste('estimatedcoef','coefsimVaccine.n',n,array_id,'.txt',sep="")
write.table(estcoef, file=fname6, sep="\t", row.names=F, col.names=T)
fname7 <- paste('covs','coefSimVaccine.n',n, array_id,'.txt',sep="")
write.table(covs, file=fname7, sep="\t",  col.names=T)

plot(int)
plot(int1)
plot(slope)

plot(intmarginalizeclosed)
plot(slopemarginalizeclosed)

# marginal plot

plot(s,mean(intmarginalizefull[burnin:SIM])+s*mean(slopemarginalizefull[burnin:SIM]),xlab="S(1)=s",type="l",ylab="E(T(1)-T(0)|s)",
     main=paste("CEP Plot\nEndpoint as Difference from Baseline\nMarginal, CI =",condindfit),ylim=c(-1,7));abline(v=0,h=0)

lines(s,mean(intmarginalizefull[burnin:SIM])+s*mean(slopemarginalizefull[burnin:SIM])+1.96*sqrt(sd(slopemarginalizefull)^2+s^2*
                                                                          sd(intmarginalizefull)^2+2*s*cov(slopemarginalizefull[burnin:SIM],intmarginalizefull[burnin:SIM])),lty=2, lwd=1,col="dark gray")
lines(s,mean(intmarginalizefull[burnin:SIM])+s*mean(slopemarginalizefull[burnin:SIM])-1.96*sqrt(sd(slopemarginalizefull[burnin:SIM])^2+s^2*
                                                                          sd(intmarginalizefull[burnin:SIM])^2+2*s*cov(slopemarginalizefull[burnin:SIM],intmarginalizefull[burnin:SIM])),lty=2, lwd=1,col="dark gray")

# conditional on age
plot(s,mean(intmarginalizeclosed[burnin:SIM])+s*mean(slopemarginalizeclosed[burnin:SIM]),xlab="S(1)=s",type="l",ylab="E(T(1)-T(0)|A,s)",main=
       paste("CEP Plot\nEndpoint as Difference from Baseline\nConditional on Age, CI =",condindfit),ylim=c(-5,5));abline(v=0,h=0)
lines(s,mean(intmarginalizeclosed[burnin:SIM])+s*mean(slopemarginalizeclosed[burnin:SIM])+1.96*sqrt(sd(slopemarginalizeclosed[burnin:SIM])^2+s^2*
                                                                              sd(intmarginalizeclosed[burnin:SIM])^2+2*s*cov(slopemarginalizeclosed[burnin:SIM],intmarginalizeclosed[burnin:SIM])),lty=2, lwd=1,col="dark gray")
lines(s,mean(intmarginalizeclosed[burnin:SIM])+s*mean(slopemarginalizeclosed[burnin:SIM])-1.96*sqrt(sd(slopemarginalizeclosed[burnin:SIM])^2+s^2*
                                                                              sd(intmarginalizeclosed[burnin:SIM])^2+2*s*cov(slopemarginalizeclosed[burnin:SIM],intmarginalizeclosed[burnin:SIM])),lty=2, lwd=1,col="dark gray")

plot(s,mean(intmarginalize[burnin:SIM])+s*mean(slopemarginalizeclosed[burnin:SIM]),xlab="S(1)=s",type="l",ylab="E(T(1)-T(0)|A,s)",main=
       paste("CEP Plot\nEndpoint as Difference from Baseline\nConditional on Age, CI =",condindfit),ylim=c(-5,5));abline(v=0,h=0)
lines(s,mean(intmarginalizeclosed[burnin:SIM])+s*mean(slopemarginalizeclosed[burnin:SIM])+1.96*sqrt(sd(slopemarginalizeclosed[burnin:SIM])^2+s^2*
                                                                                                      sd(intmarginalizeclosed[burnin:SIM])^2+2*s*cov(slopemarginalizeclosed[burnin:SIM],intmarginalizeclosed[burnin:SIM])),lty=2, lwd=1,col="dark gray")
lines(s,mean(intmarginalizeclosed[burnin:SIM])+s*mean(slopemarginalizeclosed[burnin:SIM])-1.96*sqrt(sd(slopemarginalizeclosed[burnin:SIM])^2+s^2*
                                                                                                      sd(intmarginalizeclosed[burnin:SIM])^2+2*s*cov(slopemarginalizeclosed[burnin:SIM],intmarginalizeclosed[burnin:SIM])),lty=2, lwd=1,col="dark gray")

# condition on both
plot(s,mean(int[burnin:SIM])+s*mean(slope[burnin:SIM]),xlab="S(1)=s",type="l",ylab="E(T(1)-T(0)|X,A,s)",main=
       paste("CEP Plot\nEndpoint as Difference from Baseline\nConditional on Age and Baseline, CI =",condindfit),ylim=c(-10,10));abline(v=0,h=0)
lines(s,mean(int[burnin:SIM])+s*mean(slope[burnin:SIM])+1.96*sqrt(sd(mean(slope[burnin:SIM])^2+s^2*sd(int[burnin:SIM])^2+2*s*cov(slope[burnin:SIM],int[burnin:SIM]))),lty=2, lwd=1,col="dark gray")
lines(s,mean(int[burnin:SIM])+s*mean(slope[burnin:SIM])-1.96*sqrt(sd(mean(slope[burnin:SIM])^2+s^2*sd(int[burnin:SIM])^2+2*s*cov(slope[burnin:SIM],int[burnin:SIM]))),lty=2, lwd=1,col="dark gray")

# based on Credible Interval
plot(s,mean(int[burnin:SIM])+s*mean(slope[burnin:SIM]),xlab="S(1)=s",type="l",ylab="E(T(1)-T(0)|X,A,s)",main=
       paste("CEP Plot\nEndpoint as Difference from Baseline\nConditional on Age and Baseline, CI =",condindfit),ylim=c(-10,10));abline(v=0,h=0)
lines(s, quantile(int[burnin:sim-1],  probs = 0.025, na.rm=T) + 
        quantile(slope[burnin:sim-1],  probs = 0.025, na.rm=T)*s
      ,lty=2, lwd=1,col="dark gray")
lines(s, quantile(int[burnin:sim-1],  probs = 0.975, na.rm=T) + 
        quantile(slope[burnin:sim-1],  probs = 0.975, na.rm=T)*s
      ,lty=2, lwd=1,col="dark gray")

# stratify
plot(s,mean(int[burnin:SIM]) + 5*mean(holdphi9[burnin:SIM] - holdphi6[burnin:SIM])+5^2*mean(holdphi92[burnin:SIM] - holdphi62[burnin:SIM]) +
              s*mean(slope[burnin:SIM]) + 22*mean(holdphi8[burnin:SIM] - holdphi5[burnin:SIM]), xlab="S(1)=s",type="l",ylab="E(T(1)-T(0)|X,A,s)",main=
       paste("CEP Plot\nEndpoint as Difference from Baseline\nConditional on Age and Baseline, CI =",condindfit),
     ylim=c(-3, 6));abline(v=0,h=0); legend("topleft", 
                                            legend = c("Age 5 Years", "Age 6 Years"), 
                                            fill = c("black", "red"),
                                            pt.cex = 2, 
                                            cex = 1.2)
lines(s, quantile(int[burnin:sim-1],  probs = 0.025, na.rm=T) + 
        quantile(slope[burnin:sim-1],  probs = 0.025, na.rm=T)*s  + 5*mean(holdphi9[burnin:SIM] - holdphi6[burnin:SIM])  + 
        5^2*mean(holdphi92[burnin:SIM] - holdphi62[burnin:SIM]) + 22*mean(holdphi8[burnin:SIM] - holdphi5[burnin:SIM])
      ,lty=2, lwd=1,col="dark gray")
lines(s, quantile(int[burnin:sim-1],  probs = 0.975, na.rm=T) + 
        quantile(slope[burnin:sim-1],  probs = 0.975, na.rm=T)*s  + 5*mean(holdphi9[burnin:SIM] - holdphi6[burnin:SIM]) + 
        5^2*mean(holdphi92[burnin:SIM] - holdphi62[burnin:SIM]) + 22*mean(holdphi8[burnin:SIM] - holdphi5[burnin:SIM]) 
      ,lty=2, lwd=1,col="dark gray")


lines(s, mean(int[burnin:SIM]) + 6*mean(holdphi9[burnin:SIM] - holdphi6[burnin:SIM])+
        6^2*mean(holdphi92[burnin:SIM] - holdphi62[burnin:SIM]) +
        s*mean(slope[burnin:SIM]) + 22*mean(holdphi8[burnin:SIM] - holdphi5[burnin:SIM]), type="l", lwd = 1, col = "red") 
lines(s, quantile(int[burnin:sim-1],  probs = 0.025, na.rm=T) + 
        6^2*mean(holdphi92[burnin:SIM] - holdphi62[burnin:SIM]) + 
        quantile(slope[burnin:sim-1],  probs = 0.025, na.rm=T)*s  + 6*mean(holdphi9[burnin:SIM] - holdphi6[burnin:SIM]) + 
        22*mean(holdphi8[burnin:SIM] - holdphi5[burnin:SIM])
      ,lty=2, lwd=1,col="pink")
lines(s, quantile(int[burnin:sim-1],  probs = 0.975, na.rm=T) + 
        6^2*mean(holdphi92[burnin:SIM] - holdphi62[burnin:SIM]) + 
        quantile(slope[burnin:sim-1],  probs = 0.975, na.rm=T)*s  + 6*mean(holdphi9[burnin:SIM] - holdphi6[burnin:SIM]) + 
        22*mean(holdphi8[burnin:SIM] - holdphi5[burnin:SIM])
      ,lty=2, lwd=1,col="pink")


# Region of X and A where intercept isn't through origin
plot(Age,mean(int[burnin:SIM]) + Age*(mean(holdphi9 - holdphi6)) + Age^2*mean(holdphi92) - Age^2*(mean(holdphi62)),
     xlab="A",type="l",ylab="E(T(1)-T(0)|A, X, s = 0)",main=
       "Region of A Space")

plot(Age,mean(int[burnin:SIM]) + Age*(mean(holdphi9 - holdphi6)) ,
     xlab="A",type="l",ylab="E(T(1)-T(0)|A, X, s = 0)",main=
       "Region of A Space")

plot(X,mean(int[burnin:SIM]) + X*(mean(holdphi8 - holdphi5)),xlab="X",type="l",ylab="E(T(1)-T(0)|X, A, s = 0)",main=
       "Region of X Space")


