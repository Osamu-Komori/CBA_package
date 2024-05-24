#' Generate boxplots about simulation analysis by Maxent, GM, rGM and Fisher
#'
#' \code{sim} returns boxplots of simulation analysis by Maxent, GM, rGM and Fisher

#' @param n.iter iteration number of simulation

#' @param rho correlation value

#' @param p number of variables

#' @param alpha0 a vector indicating the true values of parameters to be estimated

#' @param type type of random variables





#' @return boxplots of mean square erros and computational costs by Maxent, GM, rGM and Fisher

#' @examples
#' # Example usage
#' sim()



#' @export
sim <-function(n.iter=100,rho=0.5,p=50,alpha0=seq(0,0.1,length.out=p),type="poisson") {
library(ggplot2)
library(reshape2)
library(gridExtra)
set.seed(1)
n.grid=100
m=500
serr=NULL
time=NULL

for(i in 1:n.iter){

   A=est.sim(rho=rho,n.grid=n.grid,m=m,p=p,alpha0=alpha0,type=type)
  
   time=as.data.frame(rbind(time,A$time))
   serr=as.data.frame(rbind(serr,A$serr))
  
  colnames(time)=colnames(serr)=A$method


A1=ggplot(melt(serr), aes(x = variable, y = value)) +
  geom_boxplot(fill = "lightblue", color = "darkblue") +
  scale_x_discrete(labels = c("Maxent", "Gamma", "GM", expression(rGM[gamma %~~% 0]), expression(rGM[gamma==-0.5]), "Fisher")) +
  theme_minimal() +
  labs(title = "", x = "Methods", y = "squared errors") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 15))+scale_y_log10()+ggtitle(paste0("type=",type,", rho=",rho, ", i=",i))



A2=ggplot(melt(time), aes(x = variable, y = value)) +
  geom_boxplot(fill = "lightblue", color = "darkblue") +
  scale_x_discrete(labels = c("Maxent", "Gamma", "GM", expression(rGM[gamma %~~% 0]), expression(rGM[gamma==-0.5]), "Fisher")) +
  theme_minimal() +
  labs(title = "", x = "Methods", y = "relative computational costs") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 15))+scale_y_log10()



A = gridExtra::grid.arrange(A1, A2, nrow = 1, ncol = 2)
    print(A)
}




        }
