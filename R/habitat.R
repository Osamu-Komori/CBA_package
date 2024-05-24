#' Generation of estimated habitat maps
#'
#' \code{habitat} returns the estimated habitat maps by \code{\link{maxent}}, \code{\link{rgm}}, \code{\link{gm}} and \code{\link{fisher0}}.

#' @param Xb data matrix of environmental variables and bias variables.

#' @param sp index vector for presence locations.

#' @param iter maximum number of maximum iterations.


#' @return estimated habitat maps



#' @export
habitat <-function (env=Env2,sp = sp_pteridium_aquilinum
, iter = 500) 
{
     Xb=as.matrix(env[,-c(1,2)])

s.name=gsub("sp_","",deparse(substitute(sp)))

    library(ggplot2)
    library(gridExtra)
    Xb = as.matrix(Env2[, -c(1, 2)])
    make.prob = function(val) {
        A = exp(val)
        prob = A/sum(A)
        uni = unique(sort(prob))
        uni = quantile(uni, probs = 0:100/100)
        uni = uni[!is.na(uni)]
        B = rep(0, length(prob))
        not.na = !is.na(prob)
        for (i in 1:(length(uni) - 1)) {
            B[uni[i] < prob & prob <= uni[i + 1]] = sum(prob[not.na & 
                prob <= uni[i + 1] & prob >= 0])
        }
        return(B)
    }
    N = dim(Xb)[1]
    y = is.element(1:N, sp)
    A = maxent(env, sp, iter)
    beta.max = A$beta0
    alpha.max = A$alpha
    val.max = Xb %*% alpha.max
    AUC.max = auc(val.max, y)
    n.vari.max = sum(alpha.max != 0)
    iter.max = A$iter
    Env2$prob = make.prob(val.max)
    A1 = ggplot() + geom_raster(data = Env2, aes(x = lon, y = lat, 
        fill = prob)) + scale_fill_gradientn(colours = c(
        terrain.colors(50)[25:50]), limits = c(0, 1)) +ggtitle("Maxent")

    A = gm(env, sp, iter)
    beta.linear = A$beta0
    alpha.linear = A$alpha
    val.linear = Xb %*% alpha.linear
    AUC.linear = auc(val.linear, y)
    n.vari.linear = sum(alpha.linear != 0)
    iter.linear = A$iter
    Env2$prob = make.prob(val.linear)
    A2 = ggplot() + geom_raster(data = Env2, aes(x = lon, y = lat, 
        fill = prob)) + scale_fill_gradientn(colours = c(
        terrain.colors(50)[25:50]),  limits = c(0, 1)) +ggtitle("GM")

    A = rgm(env, sp, iter)
    beta.quadratic = A$beta0
    alpha.quadratic = A$alpha
    val.quadratic = Xb %*% alpha.quadratic
    AUC.quadratic = auc(val.quadratic, y)
    n.vari.quadratic = sum(alpha.quadratic != 0)
    iter.quadratic = A$iter
    Env2$prob = make.prob(val.quadratic)

    A3 = ggplot() + geom_raster(data = Env2, aes(x = lon, y = lat, 
        fill = prob)) + scale_fill_gradientn(colours = c(
        terrain.colors(50)[25:50]), limits = c(0, 1)) +ggtitle("rGM")



    A = fisher0(env, sp, iter)
    beta.fisher = A$beta0
    alpha.fisher = A$alpha
    val.fisher = Xb %*% alpha.fisher
    AUC.fisher = auc(val.fisher, y)
    n.vari.fisher = sum(alpha.fisher != 0)
    iter.fisher = A$iter
    Env2$prob = make.prob(val.fisher)
    A4 = ggplot() + geom_raster(data = Env2, aes(x = lon, y = lat, 
        fill = prob)) + scale_fill_gradientn(colours = c(
        terrain.colors(50)[25:50]), limits = c(0, 1)) +ggtitle("Fisher")
   A = gridExtra::grid.arrange(A1, A3, A2, A4, nrow = 2, ncol = 2)
   print(A)
}
