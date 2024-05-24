#' NCEAS data analysis by Maxent, GM, rGM and Fisher
#'
#' \code{nceas} returns the results of NCEAS data analysis by Maxent, GM, rGM and Fisher

#' @param region name of region in NCEAS data ("AWT","CAN","NSW","NZ","SA","SWI"). The default is "AWT"







#' @return boxplots of test AUC and computational costs by Maxent, GM, rGM and Fisher


#' @examples
#' # Example usage
#' nceas()



#' @export
nceas <-function (region="AWT") 
{
set.seed(1)
library(disdat)
library(ggplot2)
library(dplyr)
library(reshape2)
library(gridExtra)


M=AUCt=time=NULL
if(region=="NSW")
ID=1:54
if(region=="AWT")
ID=1:40
if(region=="CAN")
ID=1:20
if(region=="NZ")
ID=1:52
if(region=="SA")
ID=1:30
if(region=="SWI")
ID=1:30


for(id in ID){

if(id<10)
spID <-paste0(tolower(region),0,id)# "nsw14"
else
spID <-paste0(tolower(region),id)# "nsw14"



pr <- disPo(region)

bg <- disBg(region) # this is 10,000 random backgrounds



pr <- pr[pr$spid == spID, ] # subset the target species
m=dim(pr)[1]

training <- rbind(pr, bg)


if(region=="NSW"){
  if(id<=7)
    group="ba"
  else if(id<=15)
    group="db"
  else if(id<=17)
    group="nb"
  else if(id<=25)
    group="ot"
  else if(id<=33)
    group="ou"
  else if(id<=40)
    group="rt"
  else if(id<=46)
    group="ru"
  else
    group="sr" 
}
if(region=="AWT"){
  if(id<=20)
    group="bird"
  else
    group="plant"
}
if(region=="NSW"|region=="AWT"){
testing_env <- disEnv(region, group = group)
testing_pa <- disPa(region, group = group)
}
else{
testing_env <- disEnv(region)
testing_pa <- disPa(region)
}


if(region=="NSW")
covars <- c("cti", "disturb", "mi", "rainann", "raindq",
"rugged", "soildepth", "soilfert", "solrad",
"tempann", "topo", "vegsys")
if(region=="AWT")
covars=c("bc04", "bc05", "bc06", "bc12", "bc15", "slope", "topo", "tri")
if(region=="CAN")
covars=c("alt", "asp2", "ontprec", "ontslp", "onttemp", "ontveg", "watdist")
if(region=="NZ")
covars=c("age", "deficit", "hillshade", "mas", "mat", "r2pet", "slope", "sseas", "toxicats", "tseas", "vpd")
if(region=="SA")
covars=c("sabio12", "sabio15", "sabio17", "sabio18", "sabio2", "sabio4", "sabio5", "sabio6")
if(region=="SWI")
covars=c("bcc", "calc", "ccc", "ddeg", "nutri", "pday", "precyy", "sfroyy", "slope", "sradyy", "swb", "topo")

training <- training[, c("occ", covars)]
testing_env <- testing_env[, covars]





if(region=="NSW"){

training$vegsys <- as.factor(training$vegsys)

testing_env$vegsys <- as.factor(testing_env$vegsys)
  V=covars[covars!="vegsys"]

}
else if(region=="CAN"){

training$vegsys <- as.factor(training$ontveg)

testing_env$vegsys <- as.factor(testing_env$ontveg)
  V=covars[covars!="ontveg"]
}
else if(region=="NZ"){
training$vegsys <- as.factor(training$toxicats)

testing_env$vegsys <- as.factor(testing_env$toxicats)
  V=covars[covars!="toxicats"]
}
else if(region=="SWI"){
training$vegsys <- as.factor(training$calc)

testing_env$vegsys <- as.factor(testing_env$calc)
  V=covars[covars!="calc"]
}
else
  V=covars
for(v in V){
meanv <- mean(training[,v])
sdv <- sd(training[,v])
training[,v] <- (training[,v] - meanv) / sdv

testing_env[,v] <- (testing_env[,v] - meanv) / sdv
}


n.tr=dim(training)[1]
n.te=dim(testing_env)[1]




occ.tr=training$occ
total=as.data.frame(model.matrix(as.formula(~.-1),data=rbind(training[,-1],testing_env)))

train=total[1:n.tr,]
test=total[(n.tr+1):(n.tr+n.te),]




occ = testing_pa[,spID]
test=cbind(occ,test)

train=cbind(occ.tr,train)



  A=est.nceas(Xb=train,Xb.te=test)
  M=rbind(M,A)

  AUCt=rbind(AUCt,A$AUCt)
  time=rbind(time,A$time)
colnames(AUCt)=A$method
colnames(time)=A$method



df <- melt(as.data.frame(AUCt))


max_median <- max(tapply(df$value, df$variable, median))


A1 <- ggplot(df, aes(x = variable, y = value)) +
  geom_boxplot(fill = "lightblue", color = "darkblue") +
  geom_hline(yintercept = max_median, color = "red", linetype = "dashed", size = 1.5) +
  scale_x_discrete(labels = c("Maxent", "Gamma", "GM", expression(rGM[gamma %~~% 0]), expression(rGM[gamma == -0.5]), "Fisher")) +
  theme_minimal() +
  labs(title = "", x = "Methods", y = "Test AUC") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15))+ggtitle(paste0(region,", #species=",id))




df <- melt(as.data.frame(time))


max_median <- max(tapply(df$value, df$variable, median))


A2 <- ggplot(df, aes(x = variable, y = value)) +
  geom_boxplot(fill = "lightblue", color = "darkblue") +
  #geom_hline(yintercept = max_median, color = "red", linetype = "dashed", size = 1.5) +
  scale_x_discrete(labels = c("Maxent", "Gamma", "GM", expression(rGM[gamma %~~% 0]), expression(rGM[gamma == -0.5]), "Fisher")) +
  theme_minimal() +
  labs(title = "", x = "Methods", y = "relative computational costs") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15))+scale_y_log10()


A = gridExtra::grid.arrange(A1, A2, nrow = 1, ncol = 2)
    print(A)

}






}
