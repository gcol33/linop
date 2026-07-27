.onLoad <- function(libname, pkgname) {
  register_leaf_nodes()
  register_composite_nodes()
  invisible(NULL)
}
