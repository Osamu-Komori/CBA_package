#' Estimation function of Maxent
#'
#' \code{maxent} returns the estimated coefficients and the number of iterations.
#'
#' @param Xb data matrix of environmental variables. The default is the data matrix for Pteridium aquilinum.
#' @param sp index vector for presence locations. The default is the index vector for Pteridium aquilinum.
#' @param iter maximum number of iterations.
#' @param tau penalty term. A numeric value or "default" indicating that used for the linear feature in Maxent.
#' @references Phillips, S. J., & Dudík, M. (2008). Modeling of species distributions with Maxent: new extensions and a comprehensive evaluation. Ecography, 31, 161-175.
#' @return estimated coefficients and the number of iterations.
#' @seealso \code{\link{rgm}}, \code{\link{gm}}, \code{\link{gamma0}}, \code{\link{fisher0}}

#' @examples
#' # Example usage
#' maxent()



#' @export
maxent <-function(env=Env2,sp=sp_pteridium_aquilinum,iter=500,tau="default"){
   Xb=as.matrix(env[,-c(1,2)])
   m=length(sp)
   n=dim(Xb)[2]

  if(tau=="default"){
    if(m>=100)
      tau=0.05
    else if(m>=30)
      tau=-0.15/70*(m-30)+0.2
    else if(m>=10)
      tau=-0.8/20*(m-10)+1
    else 
      tau=1
  }
  
  A=sqrt(apply(Xb[sp,],2,var)/m)
  A[A==0]=min(A[A!=0])
  tau=tau*A




  ff=function(a,b,delta){
     A=-delta*a
     B=log(1+(exp(delta)-1)*b)
     C=tau*(abs(alpha+delta)-abs(alpha))
     return(A+B+C)
  }


  alpha=rep(0,n)

  Loss=numeric(0)
 
  relative=1
  for(i in 1:iter){
    a=c(apply(Xb[sp,],2,sum)/m)
    a[a==0]=10^{-10}
    a[a==1]=1-10^{-10}
    bb=c(exp(Xb%*%alpha)/sum(exp(Xb%*%alpha)))
    b=c(apply(Xb*bb,2,sum))

    delta1=-alpha
    delta2=delta3=rep(0,n)
    AA=(a-tau)*(1-b)/((1-a+tau)*b)>0
    delta2[AA]=log( ((a-tau)*(1-b)/((1-a+tau)*b))[AA] )
    AA=(a+tau)*(1-b)/((1-a-tau)*b)>0
    delta3[AA]=log( ((a+tau)*(1-b)/((1-a-tau)*b))[AA] )
    Fj=rbind(ff(a,b,delta1),ff(a,b,delta2),ff(a,b,delta3))
   
    A=which(Fj==min(Fj,na.rm=T),arr.ind=T)

    j=A[1,2]
    if(A[1,1]==1)
      delta=delta1[j]
    if(A[1,1]==2)
      delta=delta2[j]
    if(A[1,1]==3)
      delta=delta3[j]
    alpha[j]=alpha[j]+delta
    
    prob=c(exp(Xb[sp,]%*%alpha)/sum(exp(Xb%*%alpha)))
    Loss[i]=-sum(log(prob))/m+tau%*%abs(alpha)
    
    if(i>1){
      relative=(Loss[i-1]-Loss[i])/(1+abs(Loss[i]))
     if(relative<10^(-20)|Loss[i]>=Loss[i-1])
       break
    }
   
  }
 
  names(alpha)=colnames(Xb)

  return(list(alpha=alpha,iter=i))
}
