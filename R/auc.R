#' Calculation of the area under the ROC curve
#'

#' \code{auc} returns the value of the area under the ROC curve

#' @param val Values of estimated intensities of a habitat map, based on which the area under the ROC curve is calculated

#' @param y Class label taking False or True based on which the area under the ROC curve is calculated

#' @seealso \code{\link{qPPP}}

#' @export
auc <- function(val,y)
{ 
        y=y[!is.na(val)]
        val=val[!is.na(val)]
        n0 = sum(y == F)
        n1 = sum(y == T)
        v0=val[y==F]
        v1=val[y==T]

        ff=function(x)sum(x>v0)+1/2*sum(x==v0)
        v1=matrix(v1,n1,1)
        V=apply(v1,1,ff)

        AUC=sum(V)/(n0*n1)
        return(AUC)
}
