
#' Sample unique permutations from a partition space
#'
#' @param g number of generations
#' @param n total size of cluster
#' @param B number of permutations to draw
#' @return matrix of size g x B and the first row is all 1.  Each column is a unique permutation that sums to n and all are positive entries
sample_unique_perms <- function(g, n, B){

  if(n == 1) {
    return(matrix(1, ncol = B, nrow = 1))
  } else if(g == n){
    return(matrix(1, ncol = B, nrow = g))
  }

  parts <- partitions::restrictedparts(n = n-1, m = g-1,
                                       include.zero = FALSE)
  wts <- apply(parts, 2, function(x){
    if(length(unique(x)) == 1) return(1)
    cts <- RcppAlgos::permuteCount(sort(unique(x)), freqs = table(x))
    if(length(cts) > 1){
      stop()
    }
    return(cts)
  })

  if(length(wts) != ncol(parts)){
    stop("weights are too large numerically")
  }
  part_inds <- sample(1:ncol(parts), size = B,
                      replace = TRUE, prob = wts / sum(wts))
  drawn_parts <- parts[, part_inds, drop = FALSE]
  drawn_perms <- plyr::alply(.data = drawn_parts,
                             .margins = 2,
                             .fun = function(x){
                               if(length(sort(unique(x))) == 1){
                                 out <- matrix(rep(x[1], g-1), ncol = 1, nrow = g-1)
                               } else{
                                 out <- matrix(RcppAlgos::permuteSample(v = sort(unique(x)),
                                                                        freqs = table(x), n = 1),
                                               ncol = 1)

                               }
                               return(out)
                             })
  drawn_perms <- do.call('cbind', drawn_perms)
  drawn_perms <- rbind(rep(1, ncol(drawn_perms)), drawn_perms)
  return(drawn_perms)

}




#' Sample connections and form a tree
#'
#' @param gen_sizes vector of generation sizes
#' @return data frame with the following columns
#' \describe{
#' \item{gen}{generation number}
#' \item{n_in_gen}{number inside that generation}
#' \item{inf_id}{infector ID}
#' \item{id}{unique ID within cluster}
#' }
sample_connections <- function(gen_sizes){

    gen <- as.numeric(unlist(sapply(1:length(gen_sizes),
                                    function(ii) rep(ii, gen_sizes[ii]))))
  n_in_gen <- as.numeric(unlist(sapply(gen_sizes, function(x) 1:x)))
  inf_id <- sapply(1:length(gen), function(ii){
    cur_gen <- gen[ii]
    if(cur_gen == 1) return(NA)
    id <- sample(1:gen_sizes[cur_gen-1], size = 1)
    return(paste0((cur_gen-1), "-",id))
  })
  df <- data.frame(gen = gen, n_in_gen = n_in_gen, inf_id = inf_id,
                   stringsAsFactors = FALSE)
  df$id <- paste0(df$gen, "-", df$n_in_gen)
#  if("gen.1" %in% colnames(df)) browser()
  return(df)


}
