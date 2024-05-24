#' Estimation function of Gamma
#'
#' \code{gamma0} returns the estimated coefficients and a iteration number.

#' @param Xb data matrix of environmental variables. The default is the data matrix for Pteridium aquilinum.
#' @param sp index vector for presence locations. The default is the index vector for Pteridium aquilinum.
#' @param iter maximum number of iterations.
#' @param gamma gamma value in the loss function that must be negative.
#' @param tau penalty term. A numeric value or "default" indicating that used for the linear feature in Maxent.
#' @param mul multiplier applied to tau when "default" is used.
#' @references Phillips, S. J., & DudíLLk, M. (2008). Modeling of species distributions with Maxent: new extensions and a comprehensive evaluation. Ecography, 31, 161-175.
#' @return estimated coefficients and the number of iterations.
#' @seealso \code{\link{maxent}}, \code{\link{rgm}}, \code{\link{gm}}, \code{\link{fisher0}}

#' @examples
#' # Example usage
#' gamma0()

#' @export
gamma0 <-function(env=Env2,sp=sp_pteridium_aquilinum,iter=500,gamma=-0.001,tau="default",mul=1){
  Xb=as.matrix(env[,-c(1,2)])
  Xb.p=Xb[sp,]
  m=dim(Xb.p)[1]
  n=dim(Xb.p)[2]
  N=dim(Xb)[1]

  min.a=apply(Xb,2,min)
  max.b=apply(Xb,2,max)
  R=max.b-min.a

  if(tau=="default"){
    if(m>=100)
      tau=0.05
    else if(m>=30)
      tau=-0.15/70*(m-30)+0.2
    else if(m>=10)
      tau=-0.8/20*(m-10)+1
    else 
      tau=1
    tau=tau*mul
  }

 

  A=sqrt(apply(Xb.p,2,var)/m)
  A[A==0]=min(A[A!=0])
  beta=log(m)*tau*abs(gamma)*A

  gg=function(delta){
    bunshi=( exp(gamma*delta*max.b)-exp(gamma*delta*min.a) )/R*Zpf+( max.b*exp(gamma*delta*min.a)-min.a*exp(gamma*delta*max.b) )/R*Zp
        bunbo=(  ( exp((gamma+1)*delta*max.b)-exp((gamma+1)*delta*min.a) )/R*Zwf+( max.b*exp((gamma+1)*delta*min.a)-min.a*exp((gamma+1)*delta*max.b) )/R*Zw    )^(gamma/(gamma+1))

    A=log(bunshi)-log(bunbo)+beta*(abs(alpha+delta)-abs(alpha))
    return(A)
  }

  Loss=numeric(0)
  relative=1
  alpha=rep(0,n)
  Fj=matrix(0,6,n)
 
  for(i in 1:iter){
   
    Xb.p.alpha=as.numeric(Xb.p%*%alpha)
    Zp=sum(exp(gamma*Xb.p.alpha))
    Zpf=apply(exp(gamma*Xb.p.alpha)*Xb.p,2,sum)

    Xb.alpha=as.numeric(Xb%*%alpha)
    Zw=mean(exp((gamma+1)*Xb.alpha))
    Zwf=apply(exp((gamma+1)*Xb.alpha)*Xb,2,mean)
    
    delta1=-alpha

    first=gamma*(Zpf/Zp-Zwf/Zw)
    second=gamma^2*( (max.b+min.a)*Zpf/Zp -min.a*max.b-(Zpf/Zp)^2 )-gamma*(gamma+1)*( (max.b+min.a)*Zwf/Zw -min.a*max.b-(Zwf/Zw)^2 )

    first2=rep(0,length(first))
    first2[alpha!=0]=(first+beta*sign(alpha))[alpha!=0]
    first2[alpha==0&abs(first)>beta]=(first-beta*sign(first))[alpha==0&abs(first)>beta]
    delta2=-first2/second
    Fj=rbind(gg(delta1),gg(delta2))

    A=which(Fj==min(Fj,na.rm=T),arr.ind=T)
    j=A[1,2]

    if(A[1,1]==1)
      delta=delta1[j]
    if(A[1,1]==2)
      delta=delta2[j]
  
    alpha[j]=alpha[j]+delta

    Zp=sum(exp(gamma*as.numeric(Xb.p%*%alpha)))
    Zw=mean(exp((gamma+1)*as.numeric(Xb%*%alpha)))

    Loss[i]=log(Zp/(Zw)^(gamma/(gamma+1)))

    if(i>1){
     relative=(Loss[i-1]-Loss[i])/(1+abs(Loss[i]))
     if(relative<10^(-20)|Loss[i]>=Loss[i-1])
       break
     }
  }
 
  names(alpha)=colnames(Xb)

  return(list(alpha=alpha,iter=i))
}
