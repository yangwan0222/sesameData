
#' @import utils
.onAttach <- function(libname, pkgname) {
    packageStartupMessage("Loading sesameData.");
    if (has_internet()) {
        sesameDataCacheAll(showProgress = FALSE)
    }
}
