#' Generate the estimation results about simulation analysis by Maxent, GM, rGM and Fisher
#'
#' \code{est.sim} returns the results of simulation analysis by Maxent, GM, rGM and Fisher

#' @param Xb training data matrix 

#' @param Xb.te test data matrix 





#' @return boxplots of test AUC and computational costs by Maxent, GM, rGM and Fisher



#' @export
est.sim <-function (rho=0.5,n.grid=100,m=500,p=50,alpha0=seq(0,0.1,length.out=p),type="poisson"){
N=n.grid^2
library(MASS)
library(microbenchmark)
library(nortest)

pi_m=m/N

R <- matrix(rho,p,p)
diag(R)=1

if(type=="gaussian"){
mean <- rep(0,p) 
Xb <- mvrnorm(n = N, mu = mean, Sigma = R)
}

if(type=="poisson"){
lambda.poisson=3

Co=mvrnorm(n = N, mu = rep(0,p), Sigma =R )
U=apply(Co,2,pnorm)

Xb=apply(U,2,function(x)qpois(x,lambda=lambda.poisson))

}
if(type=="uniform"){

Co=mvrnorm(n = N, mu = rep(0,p), Sigma =R )
Xb=apply(Co,2,pnorm)
}
Min=apply(Xb,2,min)





        Xb1=cbind(1,Xb)	
        
        n = dim(Xb)[2]

lam=exp(Xb%*%alpha0)
     px=lam/(sum(lam))

 d <- rmultinom(n=m, size = 1, prob = px)
 



     sp=apply(d,2,function(x)which(x==1))
  

        m=length(sp)
        y=is.element(1:N,sp)
      
         
        Xb.p=Xb[sp,]


        
        fbar=apply(Xb,2,mean)
        fm=apply(Xb.p,2,mean)
        S=var(Xb)
        

maxent.loss=function(beta){
               
                A=-mean(Xb.p%*%beta)+log(sum(exp(Xb%*%beta)))
                return(A)
	}

maxent.loss.gr=function(beta){
                
                A=-apply(Xb.p,2,mean)+apply(as.numeric(exp(Xb%*%beta))*Xb,2,sum)/sum(exp(Xb%*%beta))
                return(A)
	}

gamma.loss=function(beta,gamma){

                
                A=sum(exp(gamma*Xb.p%*%beta))/(mean(exp((gamma+1)*Xb%*%beta)))^(gamma/(gamma+1))
                
                return(A)
               
                   
	}

gamma.loss.gr=function(beta,gamma){

                
                A1=gamma*apply(as.numeric(exp(gamma*Xb.p%*%beta))*Xb.p,2,sum)/(mean(exp((gamma+1)*Xb%*%beta)))^(gamma/(gamma+1))
                A2=gamma*sum(exp(gamma*Xb.p%*%beta))*apply(as.numeric(exp((gamma+1)*Xb%*%beta))*Xb,2,mean)/(mean(exp((gamma+1)*Xb%*%beta)))^((2*gamma+1)/(gamma+1))
                
                return(A1-A2)
               
                   
	}


rgm.loss=function(beta,gamma){

                
                A=sum(exp(gamma*t(t(Xb.p)-fbar)%*%beta-1/2*gamma*(gamma+1)*as.numeric(t(beta)%*%S%*%beta)   ))
                
                return(A)
               
                   
	}
rgm.loss.gr=function(beta,gamma){

                
                A=apply(as.numeric(exp(gamma*t(t(Xb.p)-fbar)%*%beta-1/2*gamma*(gamma+1)*as.numeric(t(beta)%*%S%*%beta))   )*(gamma*t(t(Xb.p)-fbar-(gamma+1)*as.numeric(S%*%beta))),2,sum)
                
                return(A)
               
                   
	}
        

gm.loss=function(beta){

               
		
                lam=exp(-Xb.p%*%beta+as.numeric(fbar%*%beta))
                
                return(sum(lam))
                
	}
gm.loss.gr=function(beta){

     lamd=-as.numeric(exp(-Xb.p%*%beta+as.numeric(fbar%*%beta)))*t(t(Xb.p)-fbar)
     return(apply(lamd,2,sum))
}





         ptm<- microbenchmark({   A=optim(par=rep(0,n),fn=maxent.loss,gr=maxent.loss.gr, method = "BFGS") },times = 1 )
       time.maxent=ptm$time

       


        par.maxent=A$par
        score=Xb%*%par.maxent
      
        serr.maxent=sum((par.maxent-alpha0)^2)
        


 

         ptm<- microbenchmark({  A=optim(par=rep(0,n),fn=gamma.loss,gr=gamma.loss.gr, method = "BFGS",gamma=-10^(-5)) },times = 1 )
       time.gamma=ptm$time

        


        par.gamma=A$par
        score=Xb%*%par.gamma
    
         serr.gamma=sum((par.gamma-alpha0)^2)
 

        ptm<- microbenchmark({ A=optim(par=rep(0,n),fn=rgm.loss,gr=rgm.loss.gr, method = "BFGS",gamma=-10^(-5))  },times = 1 )
       time.rgm=ptm$time

        
        par.rgm=A$par
        score=Xb%*%par.rgm
       
        serr.rgm=sum((par.rgm-alpha0)^2)
          

 ptm<- microbenchmark({ A=optim(par=rep(0,n),fn=rgm.loss,gr=rgm.loss.gr, method = "BFGS",gamma=-0.5)  },times = 1 )
       time.rgm_neghalf=ptm$time

        
        par.rgm_neghalf=A$par
        score=Xb%*%par.rgm_neghalf
     
        serr.rgm_neghalf=sum((par.rgm_neghalf-alpha0)^2)
       


 
        ptm<- microbenchmark({ A=optim(par=rep(0,n),fn=gm.loss,gr=gm.loss.gr, method = "BFGS")},times = 1 )
       time.gm=ptm$time
       

        par.gm=A$par
        score=Xb%*%par.gm
    
       serr.gm=sum((par.gm-alpha0)^2)

     
       




 Sinv=ginv(S)
       ptm<- microbenchmark({ par.fisher=Sinv%*%(fm-fbar)},times = 1 )
       time.fisher=ptm$time

        score=Xb%*%par.fisher
      
      serr.fisher=sum((par.fisher-alpha0)^2)


     
     


       
      time.maxent=time.maxent/time.fisher
      time.gamma=time.gamma/time.fisher
      time.gm=time.gm/time.fisher
      time.rgm=time.rgm/time.fisher
      time.rgm_neghalf=time.rgm_neghalf/time.fisher
      
      time.fisher=1
      


A=c(m,N,serr.maxent,time.maxent,par.maxent)

    
names(A)=c("m","n","serr","time",paste0("v",1:p))
A=rbind(A,c(m,N,serr.gamma,time.gamma,par.gamma))
A=rbind(A,c(m,N,serr.gm,time.gm,par.gm))
A=rbind(A,c(m,N,serr.rgm,time.rgm,par.rgm))
A=rbind(A,c(m,N,serr.rgm_neghalf,time.rgm_neghalf,par.rgm_neghalf))
A=rbind(A,c(m,N,serr.fisher,time.fisher,par.fisher))
A=data.frame(A)
method=c("Maxent","Gamma","GM","qGM","qGM_neghalf","Fisher")

A=data.frame(cbind(method,A))

return(A)




       
       

   


}
