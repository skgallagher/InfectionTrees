

#' Subset of Tuberculosis data
#'
#' @description  See references for more details.
#'  The data is collected between 2003-2009 in Maryland, USA.
#'  This is originally from \insertCite{xie2018}{InfectionTrees}

#' @format The object is a data frame of dimension 1137x7.  Each row is a detected case of TB.
#' \describe{
#' \item{sex}{Sex of individual (F/M/Unknown)}
<<<<<<< HEAD
#' \item{county}{County in Maryland of individual}
#' \item{race}{Race of individual (Asian/Black or African American/White)}
#' \item{spsmear}{Smear status of individual (positive/negative/unknown)}
#' \item{hivstatus}{HIV status of individual (positive/negative/unknown)}
#' \item{homeless}{Homelessness status of individual (yes/no)}
#' \item{groupr}{Cluster group}
#' \item{rel_time}{detection time relative to first individual in cluster}
=======
#' \item{race}{race of individual (Asian/Black or African American/White)}
#' \item{spsmear}{Smear status of individual (positive/negative/unknown)}
#' \item{hivstatus}{HIV status of individual (positive/negative/unknown)}
#' \item{homeless}{Homelessness status of individual (yes/no)}
#' \item{PCR.Cluster}{Cluster group (if any)}
#' \item{datecoll}{Date sputum was collected}
>>>>>>> fd400a1555a7a782265dd95dc3e8cb9f965adce0
#' }
#' @examples
#' head(tb_clean)
#' summary(tb_clean)
#'
#' @references
#' \insertAllCited{}
<<<<<<< HEAD
## "tb_clean"
=======
"tb_clean"
>>>>>>> fd400a1555a7a782265dd95dc3e8cb9f965adce0
