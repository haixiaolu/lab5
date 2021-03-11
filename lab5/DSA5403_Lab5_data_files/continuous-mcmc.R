#Task 3
# This function makes discrete simulations from a posterior with any number of h values
# n=nu of iterations
# You can embellish this function
simR<-function(n=10000, h=c(0.03344302,0.06165627)){
alpha<-c() # holds transition probs
alpha[1]<-1
u<-c() # holds uniform values
u[1]<-1
post<-c()# post sample
prop<-c() # vec of proposed states 1s and 2s
prop[1]=1 # initial state
post[1]=prop[1]
for(i in 2:n){ # starts at 2 because initial value given above
# proposal state 
prop[i]=sample(1:length(h),1,replace=TRUE)
# calculate alpha
# notice h[prop[i]] gives the actual value of h
alpha[i]=min(1,h[prop[i]]/h[post[i-1]])
# to calculate accepting proposal with prob alpha
# select a random value from  a uniform (0,1)
u[i]=runif(1)
if(u[i]<=alpha[i]){post[i]<-prop[i]}
else{post[i]<-post[i-1]}
}
res<-matrix(c(prop,u,alpha,post ),nc=4,nr=n,byrow=FALSE)
sim<-table(post)/n
# windows only works with a pc
# quartz with macs
dev.new(noRStudioGD = TRUE) # or quartz() 

barplot(sim)
postexact<-h/sum(h)
# The returned output is a list 
# Use obj$ to obtain whatever interests you
return(list(iter=res,sim=sim,postexact=postexact,post=post) )
}

tt = simR()
df = as.list(x = tt$post)
ggs(df)

 # Task 4
# What about different proposal distributions
# Again for the discrete case
simRQ<-function(n=1000,init=1, h=c(1,1),pr=c(1,1)/2,...){
alpha<-c() # holds transition probs
alpha[1]<-1
u<-c() # holds uniform values
u[1]<-1
post<-c()# post sample
prop<-c() # vec of proposed states 1s and 2s etc
prop[1]=init # initial state
post[1]=prop[1]
q<-function(x){pr[x]}
for(i in 2:n){ # starts at 2 because initial value given above
#make a sample from the proposal
sample(1:length(h),1,replace=TRUE,prob=pr)->prop[i]
#Calculate alpha adjusting for the proposal being non uniform
alpha[i]=min(1,h[prop[i]]*q(post[i-1])/(h[post[i-1]]*q(prop[i])))
# now choose the proposal with probability alpha
u[i]=runif(1)
if(u[i]<=alpha[i]){post[i]<-prop[i]}
else{post[i]=post[i-1]}
}
res<-matrix(c(prop,u,alpha,post ),nc=4,nr=n,byrow=FALSE,dimnames=list(1:n,c("prop","u","alpha","post")))
sim<-table(post)/n
postexact<-h/sum(h)
dev.new(noRStudioGD = TRUE)
barplot(sim,...)
tmp<-c()
ifelse(length(res[,1])>=20,tmp<-res[1:20,],tmp<-res)
return(list(iter=tmp,sim=sim,post=postexact) )
}
simRQ(1000,pr=c(1,5)/6)

simRQ(10000,h=c(3,1)  ,pr=c(1 ,2)/3)


## Now we shall look at the continuous case

# Task 5
### Using a beta proposal
### You can change the proposal to whatever you require
## a,b are the parameters of the Beta proposal
## a=b=1 is a uniform
simRC<-function(n=10,init=0.5, h=function(theta){dunif(theta)*dbinom(x=3,size=12,prob=theta)},a=1,b=1){
#dbeta(x, shape1, shape2, ncp = 0, log = FALSE)
alpha<-c() # holds transition probs
alpha[1]<-1
u<-c() # holds uniform values
u[1]<-1
post<-c()# post sample
prop<-c() # vec of proposed states 1s and 2s
prop[1]=init # initial state
post[1]=prop[1]
q=function(x){dbeta(x,a,b)}
for(i in 2:n){ # starts at 2 because initial value given above
rbeta(1,a,b)->prop[i]

alpha[i]=min(1,h(prop[i])*q(post[i-1])/(h(post[i-1])*q(prop[i])))
u[i]=runif(1)
ifelse(u[i]<=alpha[i],post[i]<-prop[i],post[i]<-post[i-1])
}
res<-matrix(c(prop,u,alpha,post ),nc=4,nr=n,byrow=FALSE,dimnames=list(1:n,c("prop","u","alpha","post")))
windows()
hist(post,freq=FALSE)

return(list(matrix=res,summary=summary(post)) )
}

simRC(n=10000)



##### LAB 4 ENDS HERE #######################

## Make a 3 state Bayes box
ddbbox3<-function(k13=4,k23=2,lik=dbinom(x=5,size=10,prob=c(0.3,0.6,0.7)),theta=c(0.3,0.6,0.7)){ # K=1...6
k21<-6*k23/k13
nuit=0
while( floor(k21)<k21 | k21>6 ){
nuit=nuit+1
ss<-sample(x=1:6,size=2,replace=TRUE)
k23<-ss[1]
k13<-ss[2]
k21<-6*k23/k13
}
library(xtable)
lik1<-lik[1] #i
lik2<-lik[2] #j
lik3<-lik[3] #l
prior3<-100
prior1<-6*prior3*lik3/(k13*lik1)
prior2<-6*prior3*lik3/(k23*lik2)
prior<-c(prior1,prior2,prior3)
prior<-prior/sum(prior)
#lik<-c(lik1,lik2)
h<-prior*lik
post=h/sum(h)
h.mat<-h%*%t(1/h)
t(apply(h.mat,c(1,2),function(x)min(x,1)))->aij
colnames(aij)<-theta
rownames(aij)<-theta
mat<-cbind(theta,prior,lik,h,post)
rownames(mat)<-1:length(lik)
Totals=c(NA,sum(prior),NA,sum(h),sum(post))
mat2=rbind(mat,Totals)
layout(matrix(c(1,2,2,4,3,3,4,3,3),nr=3,nc=3,byrow=TRUE))
rep(1,3)/3->prop
barplot(matrix(prop,nc=3,nr=1,byrow=TRUE,dimnames=list(c("Coin"),theta)),
ylim=c(0,max(prop)+0.2*max(prop)),las=1,main="Proposal\n Uniform",yaxt="n")
axis(2,at=round(0:2/6,3),labels=c("0","1/6","2/6"),
las=1,col.ticks="Red",lwd.ticks=3)
barplot(matrix(post,nc=3,nr=1,byrow=TRUE,dimnames=list(c("Coin"),theta)),
ylim=c(0,max(post)+0.2*max(post)),las=1,main="Posterior target\n post")
barplot(matrix(h,nc=3,nr=1,byrow=TRUE,dimnames=list(c("Coin"),theta)),
ylim=c(0,max(h)+0.2*max(h)),las=1,main="Proportional to target\n h")
barplot(t(aij),beside=TRUE, yaxt="n",main="aij,\n the probability of acceptance",cex.lab=5)
axis(2,at=round(0:6/6,3),labels=c("0","1/6","2/6","3/6","4/6","5/6","6/6"),
las=1,col.ticks="Blue",lwd.ticks=3)
P=aij*1/3
diag(P)=c(0,0,0)
diag(P)=1-apply(P,1, sum)
return(list(bbox=mat2,latex=xtable(mat2,digits=4),mat=mat,h=h,aij=aij,P=P,nuit=nuit))
}

ddbbox3(k13=5,k23=3,lik=dbinom(x=5,size=8,prob=c(0.4,0.5,0.8)),theta=c(0.4,0.5,0.8)) ->ans3states
ans3states$latex
xtable(ans3states$aij,dig=6)





# Make a coin-die bbox
# k = number of faces in the event set E for acceptance of proposal
# Lik = likelihood for 2 states of theta
# theta = two states
# h1 = small relative to h2
cdbbox<-function(k=1,lik,theta, h1="s"){ # K=1...6
# xtable is a library which has functions useful for latex output
library(xtable)
# rename the first and second components of the likelihood
lik1<-lik[1]
lik2<-lik[2]
# We will now make a prior that has the desired characteristics
# See pdf for proof of following
# if h1 small "s" then ... else ...
ifelse(h1=="s",pi1<-k/6*lik2/(lik1+k/6*lik2), pi1<-lik2/(lik2+k/6*lik1))
# sum of probs is 1
prior=c(pi1,1-pi1)
#lik<-c(lik1,lik2)
h<-prior*lik
# Bayes
post=h/sum(h)
# Make a matrix for the Bayes box
mat<-cbind(theta,prior,lik,h,post)
rownames(mat)<-1:length(lik)
Totals=c(NA,sum(prior),NA,sum(h),sum(post))
mat2=rbind(mat,Totals)

# Now make some plots useful in explaining the procedure
layout(matrix(c(1,2),nr=1,nc=2,byrow=TRUE))
barplot(matrix(c(0.5,0.5),nc=2,nr=1,byrow=TRUE,dimnames=list(c("Coin"),theta)),ylim=c(0,1),las=1,main="Proposal\n Uniform")
barplot(matrix(h,nc=2,nr=1,byrow=TRUE,dimnames=list(c("Coin"),theta)),ylim=c(0,max(h)+0.5*max(h)),las=1,main="Proportional to target\n h")
# Return a list of useful objects
return(list(bbox=mat2,latex=xtable(mat2,digits=6),mat=mat,h=h,h1=h1,k=k))
}

cdbbox(k=1,lik=dbinom(x=6,size=10,prob=c(0.5,0.8)),theta=c(0.5,0.8),h1="b")->ans 
ans$latex
names(ans)


## Task 2
# This function relies on getting h vectors that match E2
# Use cdbbox() to get the correct h vectors.
# Notice the ... - this corresponds to the ... in the barplot
# You can use additional arguments which pass to the barplot
coindie<-function(n=100, h=c(1/4,3/4),E2=c(5,6),init=1,...){
library(xtable)
dieset<-c()
dieset[1]<-"E1"
die<-function(n=1){
sample(1:6,size=n,replace=TRUE)
}
coin<-function(n=1){
sample(1:2,size=n,replace=TRUE)
}
face<-c()
alpha<-c() # holds acceptance probs
alpha[1]<-1
post<-c()# post sample
prop<-c() # vec of proposed states 1s and 2s
prop[1]=init # initial state
post[1]=prop[1]
dice<-c()
dice[1]<-die()

for(i in 2:n){ # starts at 2 because initial value given above
prop[i]<-coin()
alpha[i]=min(1,h[prop[i]]/h[post[i-1]])

dice[i]<-die()
ifelse(alpha[i]==1,dieset[i]<-"E1",dieset[i]<-"E2")
# is x an element of set y
if(alpha[i]==1 | (is.element(dice[i],E2) & alpha[i]!=1)){post[i]<-prop[i]}
else{post[i]<-post[i-1]}
 }  
res<-matrix(c(prop,round(alpha,2),dieset,dice,post ),nc=5,nr=n,byrow=FALSE,dimnames=list(1:n,c("proposal","alpha", "E","dice","post")))
sim<-table(post)/n
postexact<-h/sum(h)
barplot(sim,...)
return(list(iter=res,sim=sim,postexact=postexact,post=post,xtable=xtable(res,dig=1)) )
}

coindie(n=200,h=c(0.6,0.4),E2=c(3,4,5,6)) ->ans
windows()
plot(1:200, ans$post, type = "b")

coindie(n=200,h=c(1,6),E2=c(3)) ->ans
ans$it

# Use cdbbox in conjunction with coindie()
cdbbox(k=2,lik=dbinom(x=6,size=10,prob=c(0.5,0.8)),theta=c(0.5,0.8),h1="s")->ans 
names(ans)
ans2=coindie(n=30,h=ans$h,E2=1:ans$k)
ans2$it


# retrieving posterior theta values
#trace plot

ans3=coindie(n=100,h=ans$h,E2=1:ans$k)
ans4=c(0.5,0.8)[ans3$post]
plot(ans4,type="l",main="Trace plot",xlab="Iteration",ylab="theta")
