
.sesameDataGet <- function(title) {
    eh <- ExperimentHub(localHub=TRUE);
    eh <- query(eh, 'sesameData')
    if (title %in% eh$title) {
        return(eh[[which(eh$title == title)]]);
    }
    return(NULL);
}


#' Get SeSAMe data
#'
#' @param title title of the data
#' @param verbose whether to output ExperimentHub message
#' @return data object
#' @import ExperimentHub
#' @import AnnotationHub
#' @examples
#' 
#' result <- sesameDataGet('genomeInfo.hg38')
#' @export
sesameDataGet <- function(title, verbose=FALSE) {
    suppressMessages(
        log <- capture.output(
        obj <- .sesameDataGet(title)));
    obj
}

#' List all SeSAMe data
#'
#' @return all titles from SeSAMe Data
#' @examples
#' sesameDataList()
#' @export
sesameDataList <- function() {
    eh <- query(ExperimentHub(), 'sesameData')
    eh$title
}

#' Cache all SeSAMe data
#'
#' @return TRUE
#' @import ExperimentHub
#' @import AnnotationHub
#' @examples
#' sesameDataCacheAll()
#' @export
sesameDataCacheAll <- function() {
    setExperimentHubOption(arg="MAX_DOWNLOADS", 30)
    eh <- query(ExperimentHub(), 'sesameData')
    cache(eh)
    TRUE
}
