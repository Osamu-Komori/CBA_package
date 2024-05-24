#' Generate the estimation results about NCEAS data analysis by Maxent, GM, rGM and Fisher
#'
#' \code{est.nceas} returns the results of NCEAS data analysis by Maxent, GM, rGM and Fisher

#' @param Xb training data matrix 

#' @param Xb.te test data matrix 





#' @return boxplots of test AUC and computational costs by Maxent, GM, rGM and Fisher



#' @export
est.nceas <-function (Xb=nsw14.tr,Xb.te=nsw14.te){


library(MASS)
library(microbenchmark)






        Xb.p=Xb[Xb$occ==1,]
       
     
        
        n = dim(Xb)[2]-1
        m=sum(Xb$occ==1)


        N=dim(Xb)[1] 
 






     sp=which(Xb$occ==1)
    
     
        y=is.element(1:N,sp)
        y.pa=(Xb.te$occ==1)


      


Xb=as.matrix(Xb[,-1])
Xb.p=as.matrix(Xb.p[,-1])
Xb.te=as.matrix(Xb.te[,-1])

        fbar=apply(Xb,2,mean)
        fm=apply(Xb.p,2,mean)

S=var(Xb)
Sm=var(Xb.p)
vari=colnames(Xb)


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
        AUC.maxent=auc(score,y)



        score.t=Xb.te%*%par.maxent

        AUCt.maxent=auc(score.t,y.pa)
        count.maxent=A$count[2]

       
      


 

         ptm<- microbenchmark({  A=optim(par=rep(0,n),fn=gamma.loss,gr=gamma.loss.gr, method = "BFGS",gamma=-10^(-5)) },times = 1 )
       time.gamma=ptm$time

        


        par.gamma=A$par
        score=Xb%*%par.gamma
        AUC.gamma=auc(score,y)
         score.t=Xb.te%*%par.gamma
        AUCt.gamma=auc(score.t,y.pa)
         count.gamma=A$count[2]
        
        

        ptm<- microbenchmark({ A=optim(par=rep(0,n),fn=rgm.loss,gr=rgm.loss.gr, method = "BFGS",gamma=-10^(-5))  },times = 1 )
       time.rgm.loss=ptm$time

        
        par.rgm.loss=A$par
        score=Xb%*%par.rgm.loss
        AUC.rgm.loss=auc(score,y)
         score.t=Xb.te%*%par.rgm.loss
        AUCt.rgm.loss=auc(score.t,y.pa)
         count.rgm.loss=A$count[2]
       
        

 ptm<- microbenchmark({ A=optim(par=rep(0,n),fn=rgm.loss,gr=rgm.loss.gr, method = "BFGS",gamma=-0.5)  },times = 1 )
       time.rgm.loss_neghalf=ptm$time

        
        par.rgm.loss_neghalf=A$par
        score=Xb%*%par.rgm.loss_neghalf
        AUC.rgm.loss_neghalf=auc(score,y)
          score.t=Xb.te%*%par.rgm.loss_neghalf
        AUCt.rgm.loss_neghalf=auc(score.t,y.pa)
         count.rgm.loss_neghalf=A$count[2]
      
         


 
        ptm<- microbenchmark({ A=optim(par=rep(0,n),fn=gm.loss,gr=gm.loss.gr, method = "BFGS")},times = 1 )
       time.gm=ptm$time
       

        par.gm=A$par
        score=Xb%*%par.gm
        AUC.gm=auc(score,y)
         score.t=Xb.te%*%par.gm
        AUCt.gm=auc(score.t,y.pa)
        count.gm=A$count[2]
      

 
       




 Sinv=ginv(S)
       ptm<- microbenchmark({ par.fisher=Sinv%*%(fm-fbar)},times = 1 )
       time.fisher=ptm$time

        score=Xb%*%par.fisher
        AUC.fisher=auc(score,y)
         score.t=Xb.te%*%par.fisher
        AUCt.fisher=auc(score.t,y.pa)
         count.fisher=1
     

      

     
      


       
      time.maxent=time.maxent/time.fisher
      time.gamma=time.gamma/time.fisher
      time.gm=time.gm/time.fisher
      time.rgm.loss=time.rgm.loss/time.fisher
      time.rgm.loss_neghalf=time.rgm.loss_neghalf/time.fisher
     
      time.fisher=1
      


A=c(m,N,AUC.maxent,AUCt.maxent,time.maxent,count.maxent,par.maxent)

    
names(A)=c("m","n","AUC","AUCt","time","iter",vari)
A=rbind(A,c(m,N,AUC.gamma,AUCt.gamma,time.gamma,count.gamma,par.gamma))
A=rbind(A,c(m,N,AUC.gm,AUCt.gm,time.gm,count.gm,par.gm))
A=rbind(A,c(m,N,AUC.rgm.loss,AUCt.rgm.loss,time.rgm.loss,count.rgm.loss,par.rgm.loss))
A=rbind(A,c(m,N,AUC.rgm.loss_neghalf,AUCt.rgm.loss_neghalf,time.rgm.loss_neghalf,count.rgm.loss_neghalf,par.rgm.loss_neghalf))
A=rbind(A,c(m,N,AUC.fisher,AUCt.fisher,time.fisher,count.fisher,par.fisher))
A=data.frame(A)
method=c("Maxent","Gamma","GM","qGM","qGM_neghalf","Fisher")

A=data.frame(cbind(method,A))




return(A)





       
       

   


}
