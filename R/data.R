
latest_date = "2018-09-19"
cacheEnv <- new.env()

.sesameDataGet <- function(title, dateAdded = latest_date) {
    key <- paste0(title,'|',dateAdded)
    if (!exists(key, envir=cacheEnv, inherits=FALSE)) {
        eh <- query(ExperimentHub(localHub=TRUE), 'sesameData')
        obj_id <- which(eh$title == title & eh$rdatadateadded == dateAdded)
        if (length(obj_id)==1) {
            assign(key, eh[[obj_id]], envir=cacheEnv)
        } else {
            # maybe it's older version, try cache that file
            # it doesn't work properly under parallel
            # one can sesameDataCacheAll(prev_date) to avoid issue
            eh <- query(ExperimentHub(localHub=FALSE), 'sesameData')
            obj_id <- which(eh$title == title & eh$rdatadateadded == dateAdded)
            if (length(obj_id)==1) {
                cache(eh[obj_id])
                assign(key, eh[[obj_id]], envir=cacheEnv)
            } else {
                stop(
                    sprintf("%s doesn't exist. Try: sesameDataCacheAll(\"%s\")",
                    key, dateAdded))
            }
        }
    }
    return(get(key, envir=cacheEnv, inherits=FALSE))
}

#' Get SeSAMe data
#'
#' @param title title of the data
#' @param dateAdded version of the data by date added
#' @param verbose whether to output ExperimentHub message
#' @return data object
#' @import ExperimentHub
#' @import AnnotationHub
#' @examples
#' 
#' result <- sesameDataGet('genomeInfo.hg38')
#' @export
sesameDataGet <- function(title, verbose=FALSE, dateAdded = "2018-09-19") {
    if (verbose) {
        .sesameDataGet(title, dateAdded = dateAdded)
    } else {
        suppressMessages(
            log <- capture.output(
                obj <- .sesameDataGet(title, dateAdded = dateAdded)));
        obj
    }
}

#' @import curl
has_internet <- function(){
    !is.null(curl::nslookup("r-project.org", error = FALSE))
}

#' List all SeSAMe data
#'
#' @param dateAdded version of the data by date added, if "all", show all dates
#' @return all titles from SeSAMe Data
#' @examples
#' sesameDataList()
#' @export
sesameDataList <- function(dateAdded = latest_date) {
    if (has_internet()) {
        eh <- query(ExperimentHub(), 'sesameData')
    } else {
        eh <- query(ExperimentHub(localHub = TRUE), 'sesameData')
    }
    if (dateAdded == "all") {
        eh$title
    } else {
        eh$title[eh$rdatadateadded == dateAdded]
    }
}

#' List all versions of SeSAMe data
#' 
#' @return sorted unique dates of SeSAMe Data
#' @examples
#' sesameDataListDates()
#' @export
sesameDataListDates <- function() {
    if (has_internet()) {
        eh <- query(ExperimentHub(), 'sesameData')
    } else {
        eh <- query(ExperimentHub(localHub = TRUE), 'sesameData')
    }
    sort(unique(eh$rdatadateadded))
}

#' Cache all SeSAMe data
#'
#' @param dateAdded version of the data by date added, if "all", cache all dates
#' @return TRUE
#' @import ExperimentHub
#' @import AnnotationHub
#' @examples
#' sesameDataCacheAll()
#' @export
sesameDataCacheAll <- function(dateAdded = latest_date) {
    setExperimentHubOption(arg="MAX_DOWNLOADS", 100)
    eh <- query(ExperimentHub(), 'sesameData')
    if (dateAdded != "all") {
        eh <- eh[eh$rdatadateadded == dateAdded]
    }
    cache(eh)
    TRUE
}

