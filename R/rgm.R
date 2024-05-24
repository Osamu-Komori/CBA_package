#' Estimation function of rGM
#'
#' \code{rgm} returns the estimated coefficients and a iteration number.

#' @param Xb data matrix of environmental variables. The default is the data matrix for Pteridium aquilinum.
#' @param sp index vector for presence locations. The default is the index vector for Pteridium aquilinum.
#' @param iter maximum number of iterations.
#' @param gamma gamma value in the loss function that must be negative.
#' @param tau penalty term. A numeric value or "default" indicating that used for the linear feature in Maxent.
#' @param mul multiplier applied to tau when "default" is used.
#' @references Phillips, S. J., & DudíLk, M. (2008). Modeling of species distributions with Maxent: new extensions and a comprehensive evaluation. Ecography, 31, 161-175.
#' @return estimated coefficients and the number of iterations.
#' @seealso \code{\link{maxent}}, \code{\link{gm}}, \code{\link{gamma0}}, \code{\link{fisher0}}


#' @examples
#' # Example usage
#' rgm()

#' @export
rgm <-function(env=Env2,sp=sp_pteridium_aquilinum,iter=500,gamma=-0.5,tau="default",mul=0.1){
  Xb=as.matrix(env[,-c(1,2)])
   y=is.element(1:dim(Xb)[1],sp)
   Xb.p=Xb[sp,]
   m=dim(Xb.p)[1]
   n=dim(Xb.p)[2]

   fbar=apply(Xb,2,mean)
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
  beta=m*tau*abs(gamma)*A

  gg=function(Ma,Mb,L,delta){
    Cg=exp(-gamma*alphafbar-gamma*(gamma+1)/2*alphaSalpha-gamma*fbar*delta-gamma*(gamma+1)*S%*%alpha*delta-gamma*(gamma+1)/2*diag(S)*delta^2)
    A=Cg*exp(gamma*max.b*delta)/R*Ma+Cg*exp(gamma*min.a*delta)/R*Mb-L+beta*(abs(alpha+delta)-abs(alpha))
    return(t(A))
  }

  Loss=numeric(0)
  relative=1
  alpha=rep(0,n)
  Fj=matrix(0,6,n)

  for(i in 1:iter){
    alphafbar=as.numeric(alpha%*%fbar)
    alphaSalpha=as.numeric(t(alpha)%*%S%*%alpha)
    L=sum(exp(  gamma*as.numeric(t(t(Xb.p)-fbar)%*%alpha)-gamma*(gamma+1)/2*alphaSalpha  ))
    Ma=apply(exp(gamma*as.numeric(Xb.p%*%alpha))*t(t(Xb.p)-min.a),2,sum)
    Mb=apply(exp(gamma*as.numeric(Xb.p%*%alpha))*t(max.b-t(Xb.p)),2,sum)

    delta1=-alpha
  
    aj=min.a*gamma-fbar*gamma-gamma*(1+gamma)*S%*%alpha
    bj=max.b*gamma-fbar*gamma-gamma*(1+gamma)*S%*%alpha
    Aj=-1/2*gamma*(2*alphafbar+(1+gamma)*alphaSalpha)+aj^2/(2*gamma*(1+gamma)*diag(S))
    Bj=-1/2*gamma*(2*alphafbar+(1+gamma)*alphaSalpha)+bj^2/(2*gamma*(1+gamma)*diag(S))

    ggg=-gamma*alphafbar-gamma*(gamma+1)/2*alphaSalpha


    AAj=1/R*exp(Aj)*Mb
    BBj=1/R*exp(Bj)*Ma


   Rfirst=Ma*exp(ggg)*bj+Mb*exp(ggg)*aj
   Rsecond=Ma*exp(ggg)*bj^2-Ma*exp(ggg)*gamma*(gamma+1)*diag(S)+Mb*exp(ggg)*aj^2-Mb*exp(ggg)*gamma*(gamma+1)*diag(S)
   
   Rfirst2=rep(0,length(Rfirst))
   Rfirst2[alpha!=0]=(Rfirst+R*beta*sign(alpha))[alpha!=0]
   Rfirst2[alpha==0&abs(Rfirst)>R*beta]=(Rfirst-R*beta*sign(Rfirst))[alpha==0&abs(Rfirst)>R*beta]
   delta2=-Rfirst2/Rsecond

   Fj=rbind(gg(Ma,Mb,L,delta1),gg(Ma,Mb,L,delta2))
   A=which(Fj==min(Fj,na.rm=T),arr.ind=T)

   j=A[1,2]
   if(A[1,1]==1)
      delta=delta1[j]
   if(A[1,1]==2)
      delta=delta2[j]
  
   alpha[j]=alpha[j]+delta
   alphaSalpha=as.numeric(t(alpha)%*%S%*%alpha)
   Loss[i]=sum(exp(  gamma*as.numeric(t(t(Xb.p)-fbar)%*%alpha)-gamma*(gamma+1)/2*alphaSalpha  ))

   if(i>1){
     relative=(Loss[i-1]-Loss[i])/(1+abs(Loss[i]))
     if(relative<10^(-20)|Loss[i]>=Loss[i-1])
       break
     }
  }
 
  names(alpha)=colnames(Xb)

  return(list(alpha=alpha,iter=i))
}
