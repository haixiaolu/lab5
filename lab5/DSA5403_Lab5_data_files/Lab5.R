# Make a coin-die bbox
# k = number of faces in the event set E for acceptance of proposal
# Lik = likelihood for 2 states of theta
# theta = two states
# h1 = small relative to h2
cdbbox<-function(k=1,lik,theta, h1="s"){ # K=1...6
# xtable is a library which has functions useful for latex output
if(!require(xtable)) install.packages(xtable)
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
windows()
layout(matrix(c(1,2),nr=1,nc=2,byrow=TRUE))
barplot(matrix(c(0.5,0.5),nc=2,nr=1,byrow=TRUE,dimnames=list(c("Coin"),theta)),ylim=c(0,1),las=1,main="Proposal\n Uniform")
barplot(matrix(h,nc=2,nr=1,byrow=TRUE,dimnames=list(c("Coin"),theta)),ylim=c(0,max(h)+0.5*max(h)),las=1,main="Proportional to target\n h")
# Return a list of useful objects
return(list(bbox=mat2,latex=xtable(mat2,digits=6),mat=mat,h=h,h1=h1,k=k))
}

cdbbox(k=2,lik=dbinom(x=6,size=10,prob=c(0.5,0.8)),theta=c(0.5,0.8),h1="s")->ans 
ans$latex
names(ans)
ans

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

coindie(n=20,h=c(0.6,0.4),E2=c(3,4,5,6)) ->ans
ans$it

coindie(n=200,h=c(1,6),E2=c(3)) ->ans
ans$it

# Use cdbbox in conjunction with coindie()
cdbbox(k=2,lik=dbinom(x=6,size=10,prob=c(0.5,0.8)),theta=c(0.5,0.8),h1="s")->ans 
names(ans)
ans2=coindie(n=300,h=ans$h,E2=1:ans$k)
ans2$it


# retrieving posterior theta values
#trace plot

ans3=coindie(n=1000,h=ans$h,E2=1:ans$k)
ans4=c(0.5,0.8)[ans3$post]
plot(ans4,type="l",main="Trace plot",xlab="Iteration",ylab="theta")
