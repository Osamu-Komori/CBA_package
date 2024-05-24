#' Generation of estimated paths of coefficients for Maxent, GM, rGM and Fisher
#'
#' \code{make.path} returns the estimated habitat maps by \code{\link{maxent}}, \code{\link{rgm}}, \code{\link{gm}} and \code{\link{fisher0}}.

#' @param Xb data matrix of environmental variables and bias variables.

#' @param sp index vector for presence locations.

#' @param iter maximum number of maximum iterations.


#' @return estimated paths of coefficients

#' @examples
#' # Example usage
#' make.path()


#' @export
make.path <-function (env=Env2,sp=sp_pteridium_aquilinum,iter=500,Tau=c(10^{-4:2})){
library(ggplot2)
  library(gridExtra)



results=data.frame()
       for(i in 1:length(Tau)){
         if(i==1)
            A=maxent(env, sp,tau=0,iter=iter)
         else
            A=maxent(env, sp,tau=Tau[i],iter=iter)
         beta.max.default=Tau[i]
         alpha.max.default=A$alpha
         results=rbind(results,cbind(colnames(env)[-c(1,2)],alpha.max.default,Tau[i]))
         
      }
      names(results)=c("rowname","alpha","beta")
      results$alpha=as.numeric(results$alpha)
      results$beta=as.numeric(results$beta)

result_labels=results[1:length(alpha.max.default),]
A1=ggplot() +
  geom_line(data = results, aes(beta, alpha, group = rowname, color = rowname), show.legend = FALSE)+
scale_x_log10()+
geom_text(data = result_labels, aes(beta, alpha, label = rowname, color = rowname),size=3,nudge_y = 0.03, nudge_x = 1, show.legend = FALSE)+
ylab(expression(hat(alpha))) +xlab(expression(tau))+ggtitle("Maxent")



        
results=data.frame()
       for(i in 1:length(Tau)){
         if(i==1)
            A=gm(env, sp,tau=0,iter=iter)
         else
            A=gm(env, sp,tau=Tau[i],iter=iter)
         beta.max.default=Tau[i]
         alpha.max.default=A$alpha
         results=rbind(results,cbind(colnames(env)[-c(1,2)],alpha.max.default,Tau[i]))
         
      }
      names(results)=c("rowname","alpha","beta")
      results$alpha=as.numeric(results$alpha)
      results$beta=as.numeric(results$beta)

result_labels=results[1:length(alpha.max.default),]
A2=ggplot() +
  geom_line(data = results, aes(beta, alpha, group = rowname, color = rowname), show.legend = FALSE)+
scale_x_log10()+
geom_text(data = result_labels, aes(beta, alpha, label = rowname, color = rowname),size=3,nudge_y = 0.03, nudge_x = 1, show.legend = FALSE)+
ylab(expression(hat(alpha))) +xlab(expression(tau))+ggtitle("GM")


        
results=data.frame()
       for(i in 1:length(Tau)){
         if(i==1)
            A=rgm(env, sp,tau=0,iter=iter,gamma=-0.5)
         else
            A=rgm(env, sp,tau=Tau[i],iter=iter,gamma=-0.5)
         beta.max.default=Tau[i]
         alpha.max.default=A$alpha
         results=rbind(results,cbind(colnames(env)[-c(1,2)],alpha.max.default,Tau[i]))
         
      }
      names(results)=c("rowname","alpha","beta")
      results$alpha=as.numeric(results$alpha)
      results$beta=as.numeric(results$beta)

result_labels=results[1:length(alpha.max.default),]
A3=ggplot() +
  geom_line(data = results, aes(beta, alpha, group = rowname, color = rowname), show.legend = FALSE)+
scale_x_log10()+
geom_text(data = result_labels, aes(beta, alpha, label = rowname, color = rowname),size=3,nudge_y = 0.03, nudge_x = 1, show.legend = FALSE)+
ylab(expression(hat(alpha))) +xlab(expression(tau))+ggtitle("rGM")


results=data.frame()
       for(i in 1:length(Tau)){
         if(i==1)
            A=fisher0(env, sp,tau=0,iter=iter)
         else
            A=fisher0(env, sp,tau=Tau[i],iter=iter)
         beta.max.default=Tau[i]
         alpha.max.default=A$alpha
         results=rbind(results,cbind(colnames(env)[-c(1,2)],alpha.max.default,Tau[i]))
         
      }
      names(results)=c("rowname","alpha","beta")
      results$alpha=as.numeric(results$alpha)
      results$beta=as.numeric(results$beta)

result_labels=results[1:length(alpha.max.default),]
A4=ggplot() +
  geom_line(data = results, aes(beta, alpha, group = rowname, color = rowname), show.legend = FALSE)+
scale_x_log10()+
geom_text(data = result_labels, aes(beta, alpha, label = rowname, color = rowname),size=3,nudge_y = 0.03, nudge_x = 1, show.legend = FALSE)+
ylab(expression(hat(alpha))) +xlab(expression(tau))+ggtitle("Fisher")


A=gridExtra::grid.arrange(A1, A3,A2,A4, nrow = 2,ncol=2)

print(A)


}
