
#' @import utils
.onAttach <- function(libname, pkgname) {
    packageStartupMessage("Loading sesameData.");
    if (has_internet()) {
        suppressMessages(log <- capture.output(
            sesameDataCacheAll())
        );
    }
}
