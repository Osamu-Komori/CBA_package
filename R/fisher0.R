#' Estimation function of Fisher
#'
#' \code{fisher0} returns the estimated coefficients and a iteration number.

#' @param Xb data matrix of environmental variables. The default is the data matrix for Pteridium aquilinum.
#' @param sp index vector for presence locations. The default is the index vector for Pteridium aquilinum.
#' @param iter maximum number of iterations.
#' @param gamma \eqn{\gamma} value that must be negative.
#' @param tau penalty term. A numeric value or "default" indicating that used for the linear feature in Maxent.



#' @references Phillips, S. J., & Dudík, M. (2008). Modeling of species distributions with Maxent: new extensions and a comprehensive evaluation. Ecography, 31, 161-175.

#' @return estimated coefficients \eqn{\lambda} and number of iteration

#' @seealso \code{\link{maxent}}, \code{\link{rgm}}, \code{\link{gm}}, \code{\link{gamma0}}

#' @examples
#' # Example usage
#' fisher0()


#' @export
fisher0 <-function(env=Env2,sp=sp_pteridium_aquilinum,iter=500,tau="default",mul=0.1){
  Xb=as.matrix(env[,-c(1,2)])
  y=is.element(1:dim(Xb)[1],sp)
  Xb.p=Xb[sp,]
  m=dim(Xb.p)[1]
  n=dim(Xb.p)[2]

  fbar=apply(Xb,2,mean)
  fm=apply(Xb.p,2,mean)
  S=var(Xb)
   
  min.a=apply(Xb.p,2,min)
  max.b=apply(Xb.p,2,max)
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
  beta=tau*A

  gg=function(L,delta){
    
    A=1/2*diag(S)*(delta+(fbar-fm+S%*%alpha)/diag(S))^2-1/2*(fbar-fm+S%*%alpha)^2/diag(S)+as.numeric((fbar-fm)%*%alpha)+1/2*alphaSalpha-L+beta*(abs(alpha+delta)-abs(alpha))
    return(t(A))
  }
  Loss=numeric(0)
  relative=1
  alpha=rep(0,n)
  Fj=matrix(0,6,n)

  for(i in 1:iter){
    alphaSalpha=as.numeric(t(alpha)%*%S%*%alpha)
    L=as.numeric((fbar-fm)%*%alpha)+1/2*alphaSalpha+beta*sum(abs(alpha))
   
    delta1=-alpha
    delta2=-(fbar-fm+S%*%alpha+beta)/diag(S)
    delta2[alpha+delta2<0]=0
    delta3=-(fbar-fm+S%*%alpha-beta)/diag(S)
    delta3[alpha+delta3>0]=0

    Fj=rbind(gg(L,delta1),gg(L,delta2),gg(L,delta3))

    A=which(Fj==min(Fj,na.rm=T),arr.ind=T)

    j=A[1,2]

    if(A[1,1]==1)
      delta=delta1[j]
    if(A[1,1]==2)
      delta=delta2[j]
    if(A[1,1]==3)
      delta=delta3[j]
    
    alpha[j]=alpha[j]+delta

    alphaSalpha=as.numeric(t(alpha)%*%S%*%alpha)

    Loss[i]=as.numeric((fbar-fm)%*%alpha)+1/2*alphaSalpha

    if(i>1){
      relative=(Loss[i-1]-Loss[i])/(1+abs(Loss[i]))
     if(relative<10^(-20)|Loss[i]>=Loss[i-1])
       break
    }
  }
 
  names(alpha)=colnames(Xb)

  return(list(alpha=alpha,iter=i))
}
