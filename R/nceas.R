#' NCEAS data analysis by Maxent, GM, rGM and Fisher
#'
#' \code{nceas} returns the results of NCEAS data analysis by Maxent, GM, rGM and Fisher

#' @param region name of region in NCEAS data ("AWT","CAN","NSW","NZ","SA","SWI"). The default is "AWT"







#' @return boxplots of test AUC and computational costs by Maxent, GM, rGM and Fisher


#' @examples
#' # Example usage
#' nceas()



#' @export
nceas <-function (region="NZ") 
{
set.seed(1)
library(disdat)
library(ggplot2)
library(dplyr)
library(reshape2)
#Region=c("AWT","CAN","NSW","NZ","SA","SWI")

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

windows()
par(mfcol=c(1,2),mar=c(10, 4, 4, 2))
for(id in ID){
print(list(id=id))
# loading the data package

# specifying the species id
if(id<10)
spID <-paste0(tolower(region),0,id)# "nsw14"
else
spID <-paste0(tolower(region),id)# "nsw14"


# loading the presence-background data
pr <- disPo(region)




bg <- disBg(region) # this is 10,000 random backgrounds



pr <- pr[pr$spid == spID, ] # subset the target species
m=dim(pr)[1]
# combine the presence and background points
training <- rbind(pr, bg)




# loading the presence-absence testing data
# species 'nsw14' is in the 'db' group, see Elith et al. (2020) for details

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

#return(head(testing_pa))

#If region is "NSW", one of "ba"7, "db"8, "nb"2, "ot"8, "ou"8, "rt"7, "ru"6, "sr"8. 
#region is "AWT" "bird"20, "plant"20. The other regions each have only one group, so group should not be specified



# uncorrelated variables - see appendix Table 1 for variables used in each region

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
#return(head(training))

# subset uncorrelated covariates
training <- training[, c("occ", covars)]
testing_env <- testing_env[, covars]



# normalize the covariates (exept vegsys which is categorical)
# *notice: not all the models are fitted on normalized data in
# the main analysis! Please check the main text.

if(region=="NSW"){
# convert the categoricals to factor
training$vegsys <- as.factor(training$vegsys)

testing_env$vegsys <- as.factor(testing_env$vegsys)
  V=covars[covars!="vegsys"]

}
else if(region=="CAN"){
# convert the categoricals to factor
training$ontveg <- as.factor(training$ontveg)

testing_env$ontveg <- as.factor(testing_env$ontveg)
  V=covars[covars!="ontveg"]
}
else if(region=="NZ"){
training$toxicats <- as.factor(training$toxicats)

testing_env$toxicats <- as.factor(testing_env$toxicats)
  V=covars[covars!="toxicats"]
}
else if(region=="SWI"){
training$calc <- as.factor(training$calc)

testing_env$calc <- as.factor(testing_env$calc)
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



  A=est.nceas(Xb=train,Xb.pa=test)
  M=rbind(M,A)

  AUCt=rbind(AUCt,A$AUCt)
  time=rbind(time,A$time)
colnames(AUCt)=A$method
colnames(time)=A$method
boxplot(AUCt,las=2)
title(paste0("test AUC: ",region,"_",id))

boxplot(time,las=2,log="y")
title(paste0("computational costs: ",region,"_",id))

}



df <- melt(as.data.frame(AUCt))

# Calculate the maximum median value across all groups
max_median <- max(tapply(df$value, df$variable, median))

# Create the plot
A1 <- ggplot(df, aes(x = variable, y = value)) +
  geom_boxplot(fill = "lightblue", color = "darkblue") +
  geom_hline(yintercept = max_median, color = "red", linetype = "dashed", size = 1.5) +
  scale_x_discrete(labels = c("Maxent", "GM", "rGM", "Fisher")) +
  theme_minimal() +
  labs(title = "", x = "Methods", y = "Test AUC") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15))



df <- melt(as.data.frame(time))

# Calculate the maximum median value across all groups
max_median <- max(tapply(df$value, df$variable, median))

# Create the plot
A2 <- ggplot(df, aes(x = variable, y = value)) +
  geom_boxplot(fill = "lightblue", color = "darkblue") +
  scale_x_discrete(labels = c("Maxent", "GM", "rGM", "Fisher")) +
  theme_minimal() +
  labs(title = "", x = "Methods", y = "relative computational costs") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15))+scale_y_log10()




A = gridExtra::grid.arrange(A1, A2, nrow = 1, ncol = 2)
    print(A)









}
