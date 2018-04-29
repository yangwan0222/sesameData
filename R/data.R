
cacheEnv <- new.env()

#' Get SeSAMe data
#'
#' @param title title of the data
#' @return data object
#' @import ExperimentHub
#' @import AnnotationHub
#' @examples
#' 
#' sesameDataGet('genomeInfo.hg38')
#' @export
sesameDataGet <- function(title) {
    if (!exists(title, envir=cacheEnv)) {
        eh <- ExperimentHub()
        eh <- query(eh, 'sesameData')
        if (title %in% eh$title) {
            assign(title, eh$ah_id[[which(eh$title == title)]], envir=cacheEnv)
        }
    }
    return(get(title, envir=cacheEnv))
}

#' Cache all SeSAMe data
#'
#' @return TRUE
#' @import ExperimentHub
#' @import AnnotationHub
#' @examples
#' 
#' sesameDataCacheAll()
#' @export
sesameDataCacheAll <- function() {
    eh <- ExperimentHub()
    eh <- query(eh, 'sesameData')
    cache(eh)
    TRUE
}
