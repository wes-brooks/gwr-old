# Simulations for adaptive lasso estimation in the autologistic model
#
# This part contains no centering or temporal components

library(lars)
source('onestepMPLE.lars')

# Generating lattice data 

m<-15				# number of rows in the lattice
n<-15				# number of columns in the lattice
beta<-c(0,3,2,1,0,0,0)		# regression parameters	
eta<-.3				# spatial autoregression parameter
p<-length(beta)			# number of covariates (including intercept)

X<-matrix(c(rep(1,m*n),rnorm(m*n*(p-1))),m*n,p)		# covariates
X<-cbind(rep(1,m*n),scale(X[,2:p]))

# Gibbs sampling to get binary responses z

A<-func.adjmat(m,n)
z<-matrix(rbinom(m*n*2000,1,1/2),m*n,2000)
for(j in 1:2000){
	for(i in 1:(m*n)){
		z[i,j]<-rbinom(1,1,exp(X[i,]%*%beta+eta*A[i,]%*%z[,j])/(1+exp(X[i,]%*%beta+eta*A[i,]%*%z[,j])))
	}
	z[,j+1]<-z[,j]
	print(j)
}


dat<-NULL
dat$m<-m
dat$n<-n
dat$beta<-beta
dat$eta<-eta
dat$X<-X
dat$A<-A
dat$z<-z[,1901:2000]

#save(dat,file=paste('data.',m,'x',n,sep=''))

par.penalized.MPLE<-matrix(0,8,100)

# Obtain estimates for the final 100 Gibbs samples

for(j in 1901:2000){
	temp<-func.onestep.mple(c(1,1,1,1,1,1,1),0,0,X,z[,j],A,seq(0,1,.01),F,0)
	par.penalized.MPLE[,j-1900]<-c(temp[[1]],temp[[2]],temp[[3]])
	print(j)
}

#save(par.penalized.MPLE,file=paste('simresults.',m,'x',n,sep=''))

####################################################################################

# Code for processing the results (no centering, no temporal component)
#load('simresults.15x15')
results.15.15<-par.penalized.MPLE

#load('data.15x15')
data.15.15<-dat


# Proportion of samples coefficients were set to zero
apply(results.15.15[1:7,]==0,1,sum)/100

# Average number of nonzero coefficients
mean(apply(results.15.15[1:7,]!=0,2,sum))

# Average number of zero coefficients
mean(apply(results.15.15[1:7,]==0,2,sum))

# Mean coefficient value
apply(results.15.15,1,mean)

# Median coefficient value
apply(results.15.15,1,median)


# Here is some code that is intended to get the standard error
# estimators from Comets and Janzura in the Journal of Applied Probability
# (1998)

z<-data.15.15$z
X<-data.15.15$X
A<-data.15.15$A
beta0<-data.15.15$beta
eta0<-data.15.15$eta
beta<-results.15.15[1:7,]
eta<-results.15.15[8,]

n<-15*15
var.est<-matrix(0,7,100)

for(k in 1:100){
	
	nonzero.index<-beta[,k]!=0
	D<-X[,nonzero.index]*c(z[,k]-(exp(X[,nonzero.index]%*%as.vector(beta[nonzero.index,k])+eta[k]*A%*%z[,k])/(1+exp(X[,nonzero.index]%*%as.vector(beta[nonzero.index,k])+eta[k]*A%*%z[,k]))))

	A.new<-A+diag(1,n)
	Jhat<-t(D)%*%A.new%*%D

	Ihat<-matrix(0,sum(nonzero.index),sum(nonzero.index))
	w<-exp(X[,nonzero.index]%*%as.vector(beta[nonzero.index,k])+eta[k]*A%*%z[,k])/(1+exp(X[,nonzero.index]%*%as.vector(beta[nonzero.index,k])+eta[k]*A%*%z[,k]))^2

	for(i in 1:n){
		Ihat<-Ihat+w[i]*(X[i,nonzero.index]%*%t(X[i,nonzero.index]))
	}

	Ihat.inv<-solve(Ihat)
	var.cov.est <- Ihat.inv%*%Jhat%*%Ihat.inv 
	var.est[nonzero.index,k]<-diag(var.cov.est)
	var.est[!nonzero.index,k]<-NA
}

# Average variance estimate (across samples)

t(t(apply(var.est,1,median,na.rm=T)[1:4]))

# Monte Carlo variances (across samples)

for(i in 1:4){
	print(var(results.15.15[i,results.15.15[i,]!=0]))
}



#######################################################################
#######################################################################

# Simulations for adaptive lasso estimation in the autologistic model
#
# This part contains centering but no temporal components


library(lars)
source('onestepMPLE.lars')

# Generating lattice data 

m<-15
n<-15
beta<-c(0,3,2,1,0,0,0)
eta<-.3
p<-length(beta)
X<-matrix(c(rep(1,m*n),rnorm(m*n*(p-1))),m*n,p)

# Maybe need 100 (or more) simulated z's

A<-func.adjmat(m,n)
mu<-exp(X%*%beta)/(1+exp(X%*%beta))
z<-matrix(0,m*n,2000)
for(j in 1:2000){
	for(i in 1:(m*n)){
		z[i,j]<-rbinom(1,1,exp(X[i,]%*%beta+eta*A[i,]%*%(z[,j]-mu))/(1+exp(X[i,]%*%beta+eta*A[i,]%*%(z[,j]-mu))))
	}
	z[,j+1]<-z[,j]
	print(j)
}

dat<-NULL
dat$m<-m
dat$n<-n
dat$beta<-beta
dat$eta<-eta
dat$X<-X
dat$A<-A
dat$z<-z[,1901:2000]
#save(dat,file=paste('data.centered.',m,'x',n,sep=''))

par.penalized.MPLE.centered<-matrix(0,8,100)
for(j in 1901:2000){
	temp<-func.onestep.mple(c(1,1,1,1,1,1,1),0,0,X,z[,j],A,seq(0,1,.01),T,0)
	par.penalized.MPLE.centered[,j-1900]<-c(temp[[1]],temp[[2]],temp[[3]])
	print(j)
}

#save(par.penalized.MPLE.centered,file=paste('simresults.centered.',m,'x',n,sep=''))


####################################################################################

# Code for processing the results (centered, but no temporal component)


#load('simresults.centered.15x15')
results.15.15<-par.penalized.MPLE.centered


#load('data.centered.15x15')
data.15.15<-dat



# Looking at 15 x 15

apply(results.15.15==0,1,sum)/100
mean(apply(results.15.15[1:7,]!=0,2,sum))
mean(apply(results.15.15[1:7,]==0,2,sum))
apply(results.15.15,1,mean)
apply(results.15.15,1,median)


z<-data.15.15$z
X<-data.15.15$X
A<-data.15.15$A
beta0<-data.15.15$beta
eta0<-data.15.15$eta
beta<-results.15.15[1:7,]
eta<-results.15.15[8,]

n<-225
var.est<-matrix(0,7,100)


# The asymptotic var-cov given by Comets and Janzura

for(k in 1:100){

	nonzero.index<-beta[,k]!=0
	u<-exp(X[,nonzero.index]%*%as.vector(beta[nonzero.index,k]))
	mu<-u/(1+u)
	dmu.dbeta<-X[,nonzero.index]*as.vector(u/(1+u)^2)
	mu.c<-exp(X[,nonzero.index]%*%as.vector(beta[nonzero.index,k])+eta[k]*A%*%(z[,k]-mu))/(1+exp(X[,nonzero.index]%*%as.vector(beta[nonzero.index,k])+eta[k]*A%*%(z[,k]-mu)))
	X.c<-t(X[,nonzero.index])-eta[k]*t(dmu.dbeta)%*%A
	v<-(z[,k]-mu.c) 
	D<-t(X.c)*c(v)

	A.new<-A+diag(1,n)
	Jhat<-t(D)%*%A.new%*%D

	Ihat<-matrix(0,sum(nonzero.index),sum(nonzero.index))


	d2mu.dbeta2<-array(0,dim=c(sum(nonzero.index),sum(nonzero.index),n))

	for(i in 1:n){
		d2mu.dbeta2[,,i]<-((u[i]-(u[i]^2))/(1+u[i])^3)*t(X[i,nonzero.index])%*%X[i,nonzero.index]
	}

	for(i in 1:n){
		Ihat<-Ihat+eta[k]*matrix.comb(d2mu.dbeta2,A[i,])*v[i]+(mu.c[i]-(mu.c[i])^2)*X.c[,i]%*%t(X.c[,i])
	}

	Ihat.inv<-solve(Ihat)
	var.cov.est <- Ihat.inv%*%Jhat%*%Ihat.inv 
	var.est[nonzero.index,k]<-diag(var.cov.est)
	var.est[!nonzero.index,k]<-NA
	print(k)
}


# Average variance estimate (across samples)

t(t(apply(var.est,1,median,na.rm=T)[1:4]))

# Monte Carlo variances (across samples)

for(i in 1:4){
	print(var(results.15.15[i,results.15.15[i,]!=0]))
}



#######################################################################
#######################################################################

# Simulations for adaptive lasso estimation in the autologistic model
#
# This part contains no centering but temporal components

# NOTE: The data looks different in the temporal model.  When we just had spatial
# 	effects, "z" was a (m*n) x 1 vector of binary responses and "X" was a 
#	(m*n) x p matrix of covariates.  With the additional temporal component, 
#	"z" becomes a (m*n) x T matrix of binary responses and "X" becomes a 
#	(m*n) x p x T three-dimensional array of covariates.

library(lars)
source('onestepMPLE.lars')


# Generating lattice data 


m<-15
n<-15
S<-1				# Number of time units of dependence in the model
T<-5				# Number of time units of data to generate
beta<-c(0,3,2,1,0,0,0)
eta<-.3
tau<-.5
p<-length(beta)
X<-matrix(c(rep(1,m*n),rnorm(m*n*(p-1))),m*n,p)
X.array<-array(X,dim=c(length(X[,1]),length(X[1,]),T))

# Maybe need 100 (or more) simulated z's

A<-func.adjmat(m,n)
z<-array(rbinom(m*n*T*2000,1,1/2),dim=c(m*n,T,2000))
for(k in 1:2000){
	for(i in 1:(m*n)){
		z[i,1,k]<-rbinom(1,1,exp(X[i,]%*%beta+eta*A[i,]%*%z[,1,k])/(1+exp(X[i,]%*%beta+eta*A[i,]%*%z[,1,k])))
	}
	for(j in 2:T){
		for(i in 1:(m*n)){
			z[i,j,k]<-rbinom(1,1,exp(X[i,]%*%beta+eta*A[i,]%*%z[,j,k]+tau*z[i,j-1,k])/(1+exp(X[i,]%*%beta+eta*A[i,]%*%z[,j,k]+tau*z[i,j-1,k])))
		}
	}
	z[,,k+1]<-z[,,k]
	print(k)
}

dat<-NULL
dat$m<-m
dat$n<-n
dat$beta<-beta
dat$eta<-eta
dat$tau<-tau
dat$X.array<-X.array
dat$A<-A
dat$z<-z[,,1901:2000]
#save(dat,file=paste('data.temporal.',m,'x',n,'x',T,sep=''))

par.penalized.MPLE<-matrix(0,9,100)
for(k in 1901:2000){
	temp<-func.onestep.mple(c(1,1,1,1,1,1,1),0,0,X.array,z[,,k],A,seq(0,1,.01),centered=F,S=1)
	par.penalized.MPLE[,k-1900]<-c(temp[[1]],temp[[2]],temp[[3]],temp[[4]])
	print(k)
}

#save(par.penalized.MPLE,file=paste('simresults.temporal.',m,'x',n,'x',T,sep=''))

####################################################################################

# Code for processing the results (no centering, but temporal component)


#load('simresults.temporal.15x15x5')
results.15.15.5<-par.penalized.MPLE

#load('data.temporal.15x15x5')
data.15.15.5<-dat




# Looking at 15 x 15 x 5

apply(results.15.15.5==0,1,sum)/100
mean(apply(results.15.15.5[1:7,]!=0,2,sum))
mean(apply(results.15.15.5[1:7,]==0,2,sum))
apply(results.15.15.5,1,mean)
apply(results.15.15.5,1,median)



# I didn't work on the variance estimates of beta for the temporal model




