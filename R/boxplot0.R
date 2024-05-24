#' Estimation of coefficients for Maxent, GM, rGM and Fisher
#'
#' \code{boxplot0} returns the estimated coefficients \eqn{\lambda} for Maxent, GM, rGM and Fisher

#' @param Xb data matrix of environmental variables and bias variables.

#' @param sp index vector for presence locations.

#' @param iter maximum number of maximum iterations.






#' @return estimated coefficients \eqn{\lambda} for Maxent, GM, rGM and Fisher

#' @examples
#' # Example usage
#' boxplot0()


#' @export
boxplot0 <-function (env=Env2,sp=sp_pteridium_aquilinum,iter = 500) 
{
    library(ggplot2)
    A = maxent(env, sp, iter)
    alpha.max = A$alpha
    A=gamma0(env, sp, iter)
    alpha.gamma = A$alpha
    A = gm(env, sp, iter)
    alpha.gm = A$alpha
    A = rgm(env, sp, iter,gamma=-10^(-5),mul=1)
    alpha.quad0 = A$alpha
    A = rgm(env, sp, iter,gamma=-0.5)
    alpha.quad05 = A$alpha
    A = fisher0(env, sp, iter)
    alpha.fisher = A$alpha
    


        

 xd=as.data.frame(rbind(alpha.max,alpha.gamma,alpha.gm,alpha.quad0,alpha.quad05,alpha.fisher))




vari=colnames(xd)

id=(apply(xd,2,sum)!=0)
vari=vari[id]

variable= rep(vari,each=6)

method= rep(c("Maxent", "Gamma", "GM", "rGM_gamma_0", "rGM_gamma_-0.5", "Fisher"),length(vari))





coefficient=unlist(xd[,vari])

data <- data.frame(variable,method,coefficient)

data$method <- factor(data$method, levels = c("Maxent", "Gamma", "GM", "rGM_gamma_0", "rGM_gamma_-0.5", "Fisher"))




A = ggplot(data, aes(fill=method, y=coefficient, x=variable)) + 
    geom_bar(position="dodge", stat="identity") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    coord_flip() +
    ylab(expression(hat(alpha)))+
    scale_fill_manual(values = c("#FF0000", "#00FF00", "#0000FF", "#FFFF00", "#FF00FF", "#00FFFF"),
                      labels = c("Maxent", "Gamma", "GM", expression(rGM[gamma %~~% 0]), expression(rGM[gamma==-0.5]), "Fisher"))


print(A)
}
