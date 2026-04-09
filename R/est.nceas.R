#' Generate the estimation results about NCEAS data analysis by Maxent, GM, rGM and Fisher
#'
#' \code{est.nceas} returns the results of NCEAS data analysis by Maxent, GM, rGM and Fisher

#' @param Xb training data matrix 

#' @param Xb.te test data matrix 





#' @return boxplots of test AUC and computational costs by Maxent, GM, rGM and Fisher



#' @export
est.nceas <-function (Xb=nsw14.tr,Xb.pa=nsw14.te,gamma.default=-10^{-5},gamma.default2=-0.5){

library(MASS)
library(microbenchmark)
library(nortest)





        Xb.p=Xb[Xb$occ==1,]
       
     
        
        n = dim(Xb)[2]-1
        m=sum(Xb$occ==1)


        N=dim(Xb)[1] 
 


pi_m=m/N



     sp=which(Xb$occ==1)
    
     
        y=is.element(1:N,sp)
        y.pa=(Xb.pa$occ==1)


      


Xb=as.matrix(Xb[,-1])
Xb.p=as.matrix(Xb.p[,-1])
Xb.pa=as.matrix(Xb.pa[,-1])

        fbar=apply(Xb,2,mean)
        fm=apply(Xb.p,2,mean)

S=var(Xb)
Sm=var(Xb.p)
vari=colnames(Xb)

        

gm=function(beta){

               
		
                lam=exp(-Xb.p%*%beta+as.numeric(fbar%*%beta))
               
                return(sum(lam))
                
	}
gm.gr=function(beta){

     lamd=-as.numeric(exp(-Xb.p%*%beta+as.numeric(fbar%*%beta)))*t(t(Xb.p)-fbar)
     return(apply(lamd,2,sum))
}

pgm=function(beta){

                
		
                lam=exp(-Xb.p%*%beta+as.numeric(fbar%*%beta)+1/2*as.numeric(t(beta)%*%S%*%beta))
                
                return(sum(lam))
                
	}
pgm.gr=function(beta){

     lamd=as.numeric(exp(-Xb.p%*%beta+as.numeric(fbar%*%beta)))*( -t(t(Xb.p)-fbar-as.numeric(S%*%beta)))
     return(apply(lamd,2,sum))
}



gamma_linear=function(beta,gamma=gamma.default){

          
		
                lam=exp(gamma*(Xb.p%*%beta+as.numeric(fbar%*%beta)))
                
                return(sum(lam))
                
	}
gamma_linear.gr=function(beta,gamma=gamma.default){

     lamd=as.numeric(gamma*(exp(-Xb.p%*%beta+as.numeric(fbar%*%beta))))*t(t(Xb.p)-fbar)
     return(apply(lamd,2,sum))
}


ppp=function(beta){
          
                A=-sum(Xb1[sp,]%*%beta)+mean(exp(Xb1%*%beta))
                return(A)
	}

ppp.gr=function(beta){
                A=-apply(Xb1[sp,],2,sum)+apply(as.numeric(exp(Xb1%*%beta))*Xb1,2,mean)
                return(A)
	}

maxent=function(beta){
               
                A=-mean(Xb.p%*%beta)+log(sum(exp(Xb%*%beta)))
                return(A)
	}

maxent.gr=function(beta){
                
                A=-apply(Xb.p,2,mean)+apply(as.numeric(exp(Xb%*%beta))*Xb,2,sum)/sum(exp(Xb%*%beta))
                return(A)
	}

gamma0=function(beta){
                lam=exp(Xb.p%*%beta)
                A=-sum(log(lam))+m*log(mean(exp(Xb%*%beta)))
                return(A)
	}
gamma0.gr=function(beta){
                A=-apply(Xb.p,2,sum)+m/mean(exp(Xb%*%beta))*apply(as.numeric(exp(Xb%*%beta))*Xb,2,mean)
                return(A)
	}

gamma=function(beta,gamma=gamma.default){

                
                A=sum(exp(gamma*Xb.p%*%beta))/(mean(exp((gamma+1)*Xb%*%beta)))^(gamma/(gamma+1))
                
                return(A)
               
                   
	}

gamma.gr=function(beta,gamma=gamma.default){

                
                A1=gamma*apply(as.numeric(exp(gamma*Xb.p%*%beta))*Xb.p,2,sum)/(mean(exp((gamma+1)*Xb%*%beta)))^(gamma/(gamma+1))
                A2=gamma*sum(exp(gamma*Xb.p%*%beta))*apply(as.numeric(exp((gamma+1)*Xb%*%beta))*Xb,2,mean)/(mean(exp((gamma+1)*Xb%*%beta)))^((2*gamma+1)/(gamma+1))
                
                return(A1-A2)
               
                   
	}


approx=function(beta,gamma=gamma.default){

                
                A=sum(exp(gamma*t(t(Xb.p)-fbar)%*%beta-1/2*gamma*(gamma+1)*as.numeric(t(beta)%*%S%*%beta)   ))
                
                return(A)
               
                   
	}
approx.gr=function(beta,gamma=gamma.default){

                
                A=apply(as.numeric(exp(gamma*t(t(Xb.p)-fbar)%*%beta-1/2*gamma*(gamma+1)*as.numeric(t(beta)%*%S%*%beta))   )*(gamma*t(t(Xb.p)-fbar-(gamma+1)*as.numeric(S%*%beta))),2,sum)
                
                return(A)
               
                   
	}

 if(gamma.default>0)
            fnscale=-1
        else
            fnscale=1



ptm<- microbenchmark({   A=optim(par=rep(0,n),fn=maxent,gr=maxent.gr, method = "BFGS") },times = 1 )

       time.maxent=ptm$time

       


        par.maxent=A$par
        score=Xb%*%par.maxent
        AUC.maxent=auc(score,y)



        score.t=Xb.pa%*%par.maxent

        AUCt.maxent=auc(score.t,y.pa)
        count.maxent=A$count[2]

       
        p.value.maxent=lillie.test(score)$p.value


 
         
        


 ptm<- microbenchmark({ A=optim(par=rep(0,n),fn=approx,gr=approx.gr, method = "BFGS",gamma=gamma.default2,control = list(fnscale=fnscale))  },times = 1 )
       time.approx_neghalf=ptm$time

        
        par.approx_neghalf=A$par
        score=Xb%*%par.approx_neghalf
        AUC.approx_neghalf=auc(score,y)
          score.t=Xb.pa%*%par.approx_neghalf
        AUCt.approx_neghalf=auc(score.t,y.pa)
         count.approx_neghalf=A$count[2]
      
         p.value.approx_neghalf=lillie.test(score)$p.value


 
        ptm<- microbenchmark({ A=optim(par=rep(0,n),fn=gm,gr=gm.gr, method = "BFGS")},times = 1 )
       time.gm=ptm$time
       

        par.gm=A$par
        score=Xb%*%par.gm
        AUC.gm=auc(score,y)
         score.t=Xb.pa%*%par.gm
        AUCt.gm=auc(score.t,y.pa)
        count.gm=A$count[2]
      

       p.value.gm=lillie.test(score)$p.value
       




 Sinv=ginv(S)
       ptm<- microbenchmark({ par.fisher=Sinv%*%(fm-fbar)},times = 1 )
       time.fisher=ptm$time

        score=Xb%*%par.fisher
        AUC.fisher=auc(score,y)
         score.t=Xb.pa%*%par.fisher
        AUCt.fisher=auc(score.t,y.pa)
         count.fisher=1
     

      p.value.fisher=lillie.test(score)$p.value

     
      ptm<- microbenchmark({ par.fisher_m=ginv(S+pi_m*Sm)%*%(fm-fbar)},times = 1 )
       time.fisher_m=ptm$time

        score=Xb%*%par.fisher_m
        AUC.fisher_m=auc(score,y)
         score.t=Xb.pa%*%par.fisher_m
        AUCt.fisher_m=auc(score.t,y.pa)
         count.fisher_m=1
     

      p.value.fisher_m=lillie.test(score)$p.value


       
      time.maxent=time.maxent/time.fisher
      
      time.gm=time.gm/time.fisher
     
      time.approx_neghalf=time.approx_neghalf/time.fisher
      time.fisher_m=time.fisher_m/time.fisher
      time.fisher=1
      


A=c(m,N,AUC.maxent,AUCt.maxent,time.maxent,count.maxent,p.value.maxent,par.maxent)

    
names(A)=c("m","n","AUC","AUCt","time","iter","p.value",vari)

A=rbind(A,c(m,N,AUC.gm,AUCt.gm,time.gm,count.gm,p.value.gm,par.gm))

A=rbind(A,c(m,N,AUC.approx_neghalf,AUCt.approx_neghalf,time.approx_neghalf,count.approx_neghalf,p.value.approx_neghalf,par.approx_neghalf))
A=rbind(A,c(m,N,AUC.fisher,AUCt.fisher,time.fisher,count.fisher,p.value.fisher,par.fisher))

A=data.frame(A)
method=c("Maxent","GM","rGM","Fisher")

A=data.frame(cbind(method,A))




return(A)





       
       

   


}
