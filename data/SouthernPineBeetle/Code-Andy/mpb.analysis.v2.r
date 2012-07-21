library(lattice)

mpb.dat<-read.csv('mpb.csv') 

# This is the lattice

plot(Y~X,mpb.dat[mpb.dat$year==1972,],type='n')
text(mpb.dat$X,mpb.dat$Y,1:469,cex=.5)

# I have a chunk of code to get some other plots.  Let me know if
# you'd like to see them.

#######################################################


# 1 observation with NA
mpb.dat[150+469*c(0:14),]

# Remove the missing point for all years (15 years total)
mpb.dat<-mpb.dat[!mpb.dat$row%in%(150+469*c(0:14)),]


# The following code just gets the distance matrix.  There is a faster way
# with matrix multiplication but it isn't too time-consuming as it is.

dist.mat<-matrix(0,length(unique(mpb.dat$Y)),length(unique(mpb.dat$Y)))

for(i in 1:(dim(dist.mat)[[1]]-1)){
  for(j in i:dim(dist.mat)[[1]]){
    dist.mat[i,j]<-sqrt((mpb.dat$X[i]-mpb.dat$X[j])^2+(mpb.dat$Y[i]-mpb.dat$Y[j])^2)
    dist.mat[j,i]<-dist.mat[i,j]
  }
}

# A1 and A2 are different spatial neighborhoods 
# 4 nearest neighbors and 8 nearest neighbors, respectively

A1<-(dist.mat<=15000)*1
A2<-(dist.mat<=19000)*1
diag(A1)<-0
diag(A2)<-0

# Covariates are columns 14 to 27 of the mpb.dat data set with the first column of
# X given as the intercept

X1<-array(0,dim=c(dim(A1)[[1]],15,length(unique(mpb.dat$year))))

yrs<-unique(mpb.dat$year)

X1[,1,]<-1

for(t in 1:(length(yrs))){
X1[,2:15,t]<-as.matrix(mpb.dat)[mpb.dat$year==yrs[t],14:27]
}

# The binary response is whether "nifestations" is greater than 0 or not

y<-matrix(0,dim(A1)[[1]],length(unique(mpb.dat$year)))

for(t in 1:(length(yrs))){
y[,t]<-mpb.dat$nifestations[mpb.dat$year==yrs[t]]>0
}

# Proportion of plots with infestation per year
apply(y,2,mean)



###########################################################


source('onestepMPLE.lars')


# Fitting original data (no centering, two time lags of dependence)



# First, I obtain the glm estimates as initial values for the adaptive
# lasso procedure

temp.X<-NULL
temp.y<-NULL
temp.sac<-NULL
temp.tac<-NULL

for(t in 3:length(y[1,])){
temp.X<-rbind(temp.X,X1[,,t])
temp.y<-c(temp.y,y[,t])
temp.sac<-c(temp.sac,A1%*%y[,t])
temp.tac<-rbind(temp.tac,y[,(t-1):(t-2)])
}

y.glm<-temp.y
X.glm<-cbind(temp.X,temp.sac,temp.tac)

glm.est<-glm(y.glm~X.glm-1,family=binomial)

par.est1<-func.onestep.mple(glm.est$coef[1:15],glm.est$coef[16],glm.est$coef[17:18],X1,y,A1,seq(0,1,.01),centered=F,S=2)


# Compare the adaptive lasso estimates to the 
# maximum pseudolikelihood estimates (MPLE)

results.matrix<-cbind(round(c(par.est1$beta0,par.est1$beta,par.est1$eta,par.est1$tau),5),round(glm.est$coef,5))

dimnames(results.matrix)[[1]]<-c('Intercept',dimnames(mpb.dat)[[2]][14:27],'spatialAR','temporalAR_1','temporalAR_2')
dimnames(results.matrix)[[2]]<-c('adaptive lasso','MPLE')

print(results.matrix)








