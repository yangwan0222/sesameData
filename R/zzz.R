
#' @import utils
.onLoad <- function(libname, pkgname) {
    suppressMessages(log <- capture.output(
        sesameDataCacheAll()));
    fl <- system.file("extdata", "metadata.csv", package=pkgname)
    titles <- utils::read.csv(fl, stringsAsFactors=FALSE)$Title
    createHubAccessors(pkgname, titles)
}
