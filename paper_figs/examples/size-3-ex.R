## SKG
## July 8, 2020
## Maximizing the likelihood of simple transmission tree (3,2)

like <- function(par){
    p_pos <- 1 / (1 + exp(- (par[1] + par[2])))
    p_neg <- 1 / (1 + exp(- par[1]))
    loglike <- log(1 - p_pos) +
        2 * log(1- p_neg) +
        log(p_pos^2 + 3 * p_neg^2 + 2*p_pos * p_neg)

    return(-loglike)
}


par <- c(0,0)


best_pars <- optim(par = par,fn = like,
                   
