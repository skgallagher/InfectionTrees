

#' Subset of Tuberculosis data
#'
#' @description  See references for more details.
#'  The data is collected between 2003-2009 in Maryland, USA.
#'  This is originally from \insertCite{xie2018}{InfectionTrees}

#' @format The object is a data frame of dimension 1137x10.  Each row is a detected case of TB.
#' \describe{
#' \item{sex}{Sex of individual (F/M/Unknown)}
#' \item{county}{county in MD}
#' \item{ageatrept}{age of individual}
#' \item{race}{race of individual (Asian/Black or African American/White)}
#' \item{spsmear}{Smear status of individual (positive/negative/unknown)}
#' \item{hivstatus}{HIV status of individual (positive/negative/unknown)}
#' \item{homeless}{Homelessness status of individual (yes/no)}
#' \item{INIT_REGIMEN_START_DATE}{date of start of treatment}
#' \item{PCR.Cluster}{Cluster group (if any)}
#' \item{datecoll}{Date sputum was collected}
#' }
#' @examples
#' head(tb_clean)
#' summary(tb_clean)
#'
#' @references
#' \insertAllCited{}
"tb_clean"
